#ifndef _THREADED_MM_READER_
#define _THREADED_MM_READER_

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include <sys/stat.h>
#include "mmio.h"
using namespace std;

#define BATCH 10000000  // 10MB
#define MAXLINELENGTH 128

template <typename IT1, typename NT1, typename IT2, typename NT2>
void push_to_vectors(vector<IT1> & rows, vector<IT1> & cols, vector<NT1> & vals, IT2 ii, IT2 jj, NT2 vv, int symmetric, bool onebased = true)
{
    if(onebased)
    {
        ii--;  /* adjust from 1-based to 0-based */
        jj--;
    }
    rows.push_back(ii);
    cols.push_back(jj);
    vals.push_back(vv);
    if(symmetric && ii != jj)
    {
        rows.push_back(jj);
        cols.push_back(ii);
        vals.push_back(vv);
    }
}

template <typename IT1, typename NT1>
void ProcessLines(vector<IT1> & rows, vector<IT1> & cols, vector<NT1> & vals, vector<string> & lines, int symmetric, int type, bool onebased = true)
{
    if(type == 0)   // real
    {
        int64_t ii, jj;
        double vv;
        for (vector<string>::iterator itr=lines.begin(); itr != lines.end(); ++itr)
        {
            // string::c_str() -> Returns a pointer to an array that contains a null-terminated sequence of characters (i.e., a C-string)
            sscanf(itr->c_str(), "%lld %lld %lg", &ii, &jj, &vv);
            push_to_vectors(rows, cols, vals, ii, jj, vv, symmetric, onebased);
        }
    }
    else if(type == 1) // integer
    {
        int64_t ii, jj, vv;
        for (vector<string>::iterator itr=lines.begin(); itr != lines.end(); ++itr)
        {
            sscanf(itr->c_str(), "%lld %lld %lld", &ii, &jj, &vv);
            push_to_vectors(rows, cols, vals, ii, jj, vv, symmetric, onebased);
        }
    }
    else if(type == 2) // pattern
    {
        int64_t ii, jj;
        for (vector<string>::iterator itr=lines.begin(); itr != lines.end(); ++itr)
        {
            sscanf(itr->c_str(), "%lld %lld", &ii, &jj);
            push_to_vectors(rows, cols, vals, ii, jj, 1, symmetric, onebased);
        }
    }
    else
    {
        cerr << "Unrecognized matrix market scalar type" << endl;
    }
    vector<string>().swap(lines);
}




void check_newline(int *bytes_read, int bytes_requested, char *buf)
{
    if ((*bytes_read) < bytes_requested) {
        // fewer bytes than expected, this means EOF
        if (buf[(*bytes_read) - 1] != '\n') {
            // doesn't terminate with a newline, add one to prevent infinite loop later
            buf[(*bytes_read) - 1] = '\n';
            cout << "Error in Matrix Market format, appending missing newline at end of file" << endl;
            (*bytes_read)++;
        }
    }
}

// updates to curpos are reflected in the caller function
bool FetchBatch(FILE * f_local, long int & curpos, long int end_fpos, bool firstcall, vector<string> & lines, int threadid)
{
    size_t bytes2fetch = BATCH;    // we might read more than needed but no problem as we won't process them
    char * buf = new char[bytes2fetch];
    char * originalbuf = buf;   // so that we can delete it later because "buf" will move
    if(firstcall)
    {
        curpos -= 1;    // first byte is to check whether we started at the beginning of a line
        bytes2fetch += 1;
    }
    int seekfail = fseek(f_local, curpos, SEEK_SET); // move the file pointer to the beginning of thread data
    if(seekfail != 0)
        cout << "fseek failed to move to " << curpos << endl;
    
    int bytes_read = fread(buf, sizeof(char), bytes2fetch, f_local);  // read byte by byte
    if(!bytes_read)
    {
        delete [] originalbuf;
        return true;    // done
    }
    check_newline(&bytes_read, bytes2fetch, buf);
    if(firstcall)
    {
        if(buf[0] == '\n')  // we got super lucky and hit the line break
        {
            buf += 1;
            bytes_read -= 1;
            curpos += 1;
        }
        else    // skip to the next line and let the preceeding thread take care of this partial line
        {
            char *c = (char*)memchr(buf, '\n', MAXLINELENGTH); //  return a pointer to the matching byte or NULL if the character does not occur
            if (c == NULL) {
                cout << "Unexpected line without a break" << endl;
            }
            int n = c - buf + 1;
            bytes_read -= n;
            buf += n;
            curpos += n;
        }
    }
    while(bytes_read > 0 && curpos < end_fpos)  // this will also finish the last line
    {
        char *c = (char*)memchr(buf, '\n', bytes_read); //  return a pointer to the matching byte or NULL if the character does not occur
        if (c == NULL) {
            delete [] originalbuf;
            return false;  // if bytes_read stops in the middle of a line, that line will be re-read next time since curpos has not been moved forward yet
        }
        int n = c - buf + 1;
        
        // string constructor from char * buffer: copies the first n characters from the array of characters pointed by s
        lines.push_back(string(buf, n-1));  // no need to copy the newline character
        bytes_read -= n;   // reduce remaining bytes
        buf += n;   // move forward the buffer
        curpos += n;
    }
    delete [] originalbuf;
    if (curpos >= end_fpos) return true;  // don't call it again, nothing left to read
    else    return false;
}

template <typename IT, typename NT>
void ThreadedMMReader(const string & filename, vector<IT> & allrows, vector<IT> & allcols, vector<NT> & allvals, IT& numrows, IT& numcols, bool & isWeighted)
{
    int32_t type = -1;
    int32_t symmetric = 0;
    int64_t nrows, ncols, nonzeros;
    int64_t linesread = 0;
    numrows = numcols = 0;
    isWeighted = true;
    
    FILE *f;
    MM_typecode matcode;
    if ((f = fopen(filename.c_str(), "r")) == NULL)
    {
        printf("Matrix-market file can not be found\n");
        exit(1);
    }
    int mm_value = mm_read_banner(f, &matcode);
    if (mm_value != 0)
    {
        if(mm_value == MM_NO_HEADER)
        {
            printf("No Matrix Market Header found, assuming Float(real), general(non-symmetric) sparse matrix\n");
            mm_set_sparse(&matcode);
            mm_set_real(&matcode);
            mm_set_general(&matcode);
            type = 0;
        }
        else
        {
            printf("Could not process Matrix Market banner.\n");
            exit(1);
        }
    }
    linesread++;
    
    if (mm_is_complex(matcode))
    {
        printf("Sorry, this application does not support complex types");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
    }
    else if(mm_is_real(matcode))
    {
        cout << "Matrix is Float" << endl;
        type = 0;
    }
    else if(mm_is_integer(matcode))
    {
        cout << "Matrix is Integer" << endl;
        type = 1;
    }
    else if(mm_is_pattern(matcode))
    {
        cout << "Matrix is Boolean" << endl;
        type = 2;
        isWeighted = false;
    }
    if(mm_is_symmetric(matcode))
    {
        cout << "Matrix is symmetric" << endl;
        symmetric = 1;
    }
    int ret_code;
    if ((ret_code = mm_read_mtx_crd_size(f, &nrows, &ncols, &nonzeros, &linesread)) !=0)  // ABAB: mm_read_mtx_crd_size made 64-bit friendly
        exit(1);
    
    numrows = static_cast<IT> (nrows);
    numcols = static_cast<IT> (ncols);
    cout << "Total number of nonzeros expected across all processors is " << nonzeros << endl;

    
    
    // Use fseek again to go backwards two bytes and check that byte with fgetc
    struct stat st;     // get file size
    if (stat(filename.c_str(), &st) == -1)
    {
        exit(1);
    }
    int64_t file_size = st.st_size;
    cout << "File is " << file_size << " bytes" << endl;
    long int ffirst = ftell(f);
    fclose(f);
    
    vector<IT> localsizes(omp_get_max_threads());
    
#pragma omp parallel
    {
        long int fpos, end_fpos;
        int this_thread = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
            
        if(this_thread == 0)
            fpos = ffirst;
        else
            fpos = this_thread * file_size / num_threads;
        
        if(this_thread != (num_threads-1)) end_fpos = (this_thread + 1) * file_size / num_threads;
        else end_fpos = file_size;
        
        
        FILE * f_perthread = fopen(filename.c_str(), "rb");   // reopen
        
        vector<string> lines;
        bool finished = FetchBatch(f_perthread, fpos, end_fpos, true, lines, this_thread);
        int64_t entriesread = lines.size();
        
        vector<IT> rows;
        vector<IT> cols;
        vector<NT> vals;
        
        ProcessLines(rows, cols, vals, lines, symmetric, type);
        while(!finished)
        {
            finished = FetchBatch(f_perthread, fpos, end_fpos, false, lines, this_thread);
            entriesread += lines.size();
            ProcessLines(rows, cols, vals, lines, symmetric, type);
        }
        localsizes[this_thread] = rows.size();
        #pragma omp barrier
        
        #pragma omp single 
        {
            size_t nnz_after_symmetry = std::accumulate(localsizes.begin(), localsizes.begin()+num_threads, IT(0));
            
            allrows.resize(nnz_after_symmetry);
            allcols.resize(nnz_after_symmetry);
            allvals.resize(nnz_after_symmetry);
            
            // copy(localsizes.begin(), localsizes.end(), ostream_iterator<IT>(cout, " ")); cout << endl;
        }
        
        IT untilnow = std::accumulate(localsizes.begin(), localsizes.begin()+this_thread, IT(0));
        
        std::copy(rows.begin(), rows.end(), allrows.begin() + untilnow);
        std::copy(cols.begin(), cols.end(), allcols.begin() + untilnow);
        std::copy(vals.begin(), vals.end(), allvals.begin() + untilnow);
    }
}

#endif


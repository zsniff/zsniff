/** File: zsniff.cpp
    Author: Candice Quates 
    Program: zsniff deflate fragment guesser/extractor 

    Takes argument of a fragment, optionally --extract outputname
    Outputs information about deflate streams detected and their
    likely contents.
*/

/* compile with:  
g++ -std=c++0x -o zsniff zdata.cpp zsniff.cpp
*/


#include <vector>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#include "zsniff.h"

#define VERSION "zsniff 0.2 by Candice Quates, Vassil Roussev 07/2013"

void loadFile(std::vector<unsigned char>& buffer, const std::string& filename) //designed for loading files from hard disk in an std::vector
{
    std::ifstream file(filename.c_str(), std::ios::in|std::ios::binary|std::ios::ate);

    //get filesize
    std::streamsize size = 0;
    if(file.seekg(0, std::ios::end).good()) size = file.tellg();
    if(file.seekg(0, std::ios::beg).good()) size -= file.tellg();

    //read contents of the file into the vector
    if(size > 0)
    {
        buffer.resize((size_t)size);
        file.read((char*)(&buffer[0]), size);
    }
    else buffer.clear();
}

int
writeExtracted(std::vector<unsigned char> &image, std::string filename)  {
    std::string outfilename = filename+".undefl";
    std::filebuf fb;
    fb.open (outfilename.c_str(),std::ios::out|std::ios::binary);
    if (fb.is_open()) {
        std::ostream os(&fb);
        os.write((const char *)image.data(), image.size());
        fb.close();
        return 0;
    } else {
        std::cerr << "zsniff: cannot write to " << outfilename << std::endl;
        return 1;
    }
}

// going back to basics, process a file, do statistics on it.
// block-based sampling mode == later
// option for extract/noextract - removed a lot of that code.
int 
main(int argc, char *argv[])
{
    bool extract = false;
    if (argc < 2 ) {
        std::cout << VERSION << std::endl << std::endl;
        std::cout << "Usage: "<< argv[0] << " [--extract] [--block 65536] filename(s) " <<  std::endl;
        return 1;
    }   
    bool blockm=false;
    int blocksize=0;
    // look for extract and block size
    for (int i=0; i<argc; i++) {
        if (std::string("--extract").compare(std::string(argv[i]))==0) {
            extract=true;
        }
        if (blockm==true) {
            if (isdigit(argv[i][0])) {
                blocksize=atoi(argv[i]); 
            } else {
                std::cerr << "usage: "<< argv[0] << " [--extract] [--block 65535] filename(s) " <<  std::endl;
                return 1;
            }
            blockm=false;
        }
        if (std::string("--block").compare(std::string(argv[i]))==0) {
            blockm=true;
        }
    }
    blockm=false;
    //loop for multiple files here 
    for (int i=1; i<argc; i++) {
        if (std::string("--extract").compare(std::string(argv[i]))==0) {
            continue;
        } 
        // block marker to ignore argument after --block
        if (blockm==true) {
            blockm=false;
            continue;
        }
        if (std::string("--block").compare(std::string(argv[i]))==0) {
            blockm=true;
            continue;
        }
        const char* filename = argv[i] ;

        //load and decode
        std::vector<unsigned char> buffer;
        loadFile(buffer, filename);
        unsigned long w, h;
        if (buffer.empty()) {
            //send nothing.
            std::cerr << "zsniff: file " << filename << " empty or not found" << std::endl;
            return 1;
        }
        int error;
        // note if blocking... 
        if (blocksize > 0 && buffer.size() > blocksize) {
            int currentoffset=0;
            while (buffer.size() >= currentoffset+blocksize) {
                std::vector<unsigned char> image;
                std::cout << "zsniff processing: "<< filename << ":"<< currentoffset << std::endl;
                error = decodeStuff(image, &buffer[currentoffset], (unsigned long)blocksize,filename);
                currentoffset+=blocksize;
                if (extract && image.size() > 0)  {
                    writeExtracted(image,std::string(filename)+std::to_string(currentoffset));
                }
            } 
            if (currentoffset < buffer.size()) { // leftovers 
                std::vector<unsigned char> image;
                std::cout << "zsniff processing: "<< filename << ":"<< currentoffset << std::endl;
                error = decodeStuff(image, &buffer[currentoffset], (unsigned long)buffer.size()-currentoffset,filename);
                if (extract && image.size() > 0)  {
                    writeExtracted(image,std::string(filename)+std::to_string(currentoffset));
                }
            }
        } else {
            std::vector<unsigned char> image;
            std::cout << "zsniff processing: "<< filename << std::endl;
            error = decodeStuff(image,&buffer[0],(unsigned long)buffer.size(),filename);

            //if there's an error, check type
            // error 52 == endcase
            if (error == 48) {
                std::cerr << "zsniff: file " << filename << " empty or not found " << std::endl;
                return 2;
            }
            if(error != 0) 
                std::cout << "zsniff: error: " << error << std::endl;
            if (extract && image.size() > 0)  {
                writeExtracted(image,std::string(filename));
            }
        }
    }
// end loop.... 
    return 0;
}




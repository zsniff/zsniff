/** File: zdata.cpp
    Maintainer: Candice Quates
    Contents: Fixed-up, slider added, otherwise mangled picoPNG code. 

    Should be standalone enough to call the decodeStuff function.
*/

#include <vector>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

// Note:  Software attacked with axe to remove all PNG functionality and
// preserve DEFLATE code.  Many bugs fixed. -- CEQ 2/1/13
// note software IS modified
// 
// picoPNG version 20101224
// Copyright (c) 2005-2010 Lode Vandevenne
//
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
//     1. The origin of this software must not be misrepresented; you must not
//     claim that you wrote the original software. If you use this software
//     in a product, an acknowledgment in the product documentation would be
//     appreciated but is not required.
//     2. Altered source versions must be plainly marked as such, and must not be
//     misrepresented as being the original software.
//     3. This notice may not be removed or altered from any source distribution.

using namespace std;

static const unsigned long LENBASE[29] =  {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258};
static const unsigned long LENEXTRA[29] = {0,0,0,0,0,0,0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,  4,  5,  5,  5,  5,  0};
static const unsigned long DISTBASE[30] =  {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};
static const unsigned long DISTEXTRA[30] = {0,0,0,0,1,1,2, 2, 3, 3, 4, 4, 5, 5,  6,  6,  7,  7,  8,  8,   9,   9,  10,  10,  11,  11,  12,   12,   13,   13};
static const unsigned long CLCL[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15}; //code length code lengths

struct HuffmanTree
{
    int makeFromLengths(const std::vector<unsigned long>& bitlen, unsigned long maxbitlen)
    { 
        //make tree given the lengths
        unsigned long numcodes = (unsigned long)(bitlen.size()), treepos = 0, nodefilled = 0;
        std::vector<unsigned long> blcount(maxbitlen + 1, 0), nextcode(maxbitlen + 1, 0);
        std::vector<unsigned long> tree1d(numcodes);
        for(unsigned long bits = 0; bits < numcodes; bits++) 
            blcount[bitlen[bits]]++; //count number of instances of each code length
        for(unsigned long bits = 1; bits <= maxbitlen; bits++) 
            nextcode[bits] = (nextcode[bits - 1] + blcount[bits - 1]) << 1;
        for(unsigned long n = 0; n < numcodes; n++) 
            if(bitlen[n] != 0) 
                tree1d[n] = nextcode[bitlen[n]]++; //generate all the codes
        tree2d.clear(); 
        tree2d.resize(numcodes * 2, 32767); //32767 here means the tree2d isn't filled there yet
        tree2d_len.clear(); 
        tree2d_len.resize(numcodes * 2, 32767); //32767 here means the tree2d isn't filled there yet
        for(unsigned long n = 0; n < numcodes; n++) //the codes
            for(unsigned long i = 0; i < bitlen[n]; i++) //the bits for this code
            {
                unsigned long bit = (tree1d[n] >> (bitlen[n] - i - 1)) & 1;
                if(treepos > numcodes - 2) return 55;
                if(tree2d[2 * treepos + bit] == 32767) //not yet filled in
                {
                    if(i + 1 == bitlen[n]) { 
                        tree2d[2 * treepos + bit] = n; 
                        tree2d_len[2 * treepos + bit] = bitlen[n]; 
                        treepos = 0; 
                    } //last bit
                    else { 
                        tree2d[2 * treepos + bit] = ++nodefilled + numcodes; 
                        treepos = nodefilled; 
                    } //addresses are encoded as values > numcodes
                }
                else {
                    treepos = tree2d[2 * treepos + bit] - numcodes; //subtract numcodes from address to get address value
                }
            }
            return 0;
    }

    int decode(bool& decoded, unsigned long& result, size_t& treepos, unsigned long bit) const
    {   //Decodes a symbol from the tree
        unsigned long numcodes = (unsigned long)tree2d.size() / 2;
        if(treepos >= numcodes) 
            return 11; //error: you appeared outside the codetree
        result = tree2d[2 * treepos + bit];
        decoded = (result < numcodes);
        treepos = decoded ? 0 : result - numcodes;
        return 0;
    }
    std::vector<unsigned long> tree2d; //2D representation of a huffman tree: The one dimension is "0" or "1", the other contains all nodes and leaves of the tree.
    std::vector<unsigned long> tree2d_len; //2D representation of a huffman tree: The one dimension is "0" or "1", the other contains all nodes and leaves of the tree.
};


struct Inflator
{
    int error;
    int blockcount;
    int xmlcount;
    int execount;
    int pngcount;
    static unsigned long readBitFromStream(size_t& bitp, const unsigned char* bits) { 
        unsigned long result = (bits[bitp >> 3] >> (bitp & 0x7)) & 1; 
        bitp++; 
        return result;
    }
    static unsigned long readBitsFromStream(size_t& bitp, const unsigned char* bits, size_t nbits)
    {
        unsigned long result = 0;
        for(size_t i = 0; i < nbits; i++) 
            result += (readBitFromStream(bitp, bits)) << i;
        return result;
    }
    void inflate(std::vector<unsigned char>& out,const std::vector<unsigned char>& in,  size_t inpos = 0)
    {
        size_t bp = 0, pos = 0; //bit pointer and byte pointer
        error = 0;
        unsigned long BFINAL = 0;
        unsigned long BTYPE = 0;
        unsigned long BTYPE2 = 0;
        int slide=0;
        size_t bpsave = 0;
        blockcount=0;
        xmlcount=0;
        execount=0;
        pngcount=0;
        int tmpcount=0;
        int tmpsize=0;
        int previous=0;
        while(!error)
        {
            // if EOF
            if(bp >> 3 >= in.size()) 
            { 
                error = 0; 
                if (blockcount==1)
                    cout << "type|block count|raw size" << endl;
                switch (previous) {
                    case 2:
                        cout << "xml|" << tmpcount << "|"<< (int)(tmpsize / 8 ) << endl;
                        break;
                    case 3:
                        cout << "exe|" << tmpcount <<"|"<<(int)(tmpsize/8) << endl;
                        break;
                    case 4:
                        cout << "png|" << tmpcount << "|"<<(int)(tmpsize/8) << endl;
                        break;
                }
                return; 
            } 
            //SLIDING CODE HERE 
            if (!slide) {
                BTYPE = readBitFromStream(bp, &in[inpos]); 
                BTYPE2 = readBitFromStream(bp, &in[inpos]); // reads type code to check.
                BTYPE += 2 * BTYPE2;
            } else {
                BTYPE = BTYPE2; // pull previous loop iteration's btype
                BTYPE2 = readBitFromStream(bp, &in[inpos]); // reads type code to check.
                BTYPE += 2 * BTYPE2;
            }// end of slide.
            if (BTYPE== 2) { 
                bpsave=bp;
                int done=inflateHuffmanBlock(out,&in[inpos], bp, pos, in.size(), BTYPE);
                //if (done) // if we are only caring about the first block we find, breaks here
                //    break;
                if (previous==0)
                    previous=done;
                switch (done) {
                   case 2:
                       xmlcount++;
                       break;
                   case 3:
                       execount++;
                       break;
                   case 4:
                       pngcount++;
                       break;
                }
                if (error) {
                    error=0; // if there was an error, reset and try again.  
                    bp=bpsave; 
                } else { 
                    if (done!=previous) {
                        if (blockcount==tmpcount)
                            cout << "type|block count|raw size" << endl;
                        switch (previous) {
                            case 2:
                                cout << "xml|" << tmpcount << "|"<<(int)(tmpsize/8) << endl;
                                break;
                            case 3:
                                cout << "exe|" << tmpcount << "|"<<(int)(tmpsize/8) << endl;
                                break;
                            case 4:
                                cout << "png|" << tmpcount << "|"<<(int)(tmpsize/8) << endl;
                                break;
                        }
                        tmpcount=0;
                        tmpsize=0;
                     }
                     blockcount++;
                     previous=done;
                     tmpcount++;
                     tmpsize+=bp-bpsave;
                 }
                    slide=0;
                } else {
                    slide=1;
                }
            }
        if(!error) 
            out.resize(pos);
    }

    HuffmanTree codetree, codetreeD, codelengthcodetree; 

    //decode a single symbol from given list of bits with given code tree. return value is the symbol
    unsigned long huffmanDecodeSymbol(const unsigned char* in, size_t& bp, const HuffmanTree& codetree, size_t inlength)
    { 
        
        bool decoded; 
        unsigned long ct;
        for(size_t treepos = 0;;)
        {
            if((bp & 0x07) == 0 && (bp >> 3) > inlength) 
            { 
                error = 10; 
                return 0; 
            } //error: end reached without endcode
            error = codetree.decode(decoded, ct, treepos, readBitFromStream(bp, in)); 
            if(error) 
                return 0; //stop, an error happened
            if(decoded) 
                return ct;
        }
    }

    //get the tree of a deflated block with dynamic tree, the tree itself is also Huffman compressed with a known tree
    void getTreeInflateDynamic(HuffmanTree& tree, HuffmanTree& treeD, const unsigned char* in, size_t& bp, size_t inlength)
    { 
        std::vector<unsigned long> bitlen(288, 0), bitlenD(32, 0);
        if(bp >> 3 >= inlength - 2) //the bit pointer is or will go past the memory
        { 
            error = 49; 
            return; 
        } 
        size_t HLIT =  readBitsFromStream(bp, in, 5) + 257; //number of literal/length codes + 257
        size_t HDIST = readBitsFromStream(bp, in, 5) + 1; //number of dist codes + 1
        size_t HCLEN = readBitsFromStream(bp, in, 4) + 4; //number of code length codes + 4
        std::vector<unsigned long> codelengthcode(19); //lengths of tree to decode the lengths of the dynamic tree
        for(size_t i = 0; i < 19; i++) 
            codelengthcode[CLCL[i]] = (i < HCLEN) ? readBitsFromStream(bp, in, 3) : 0;
        error = codelengthcodetree.makeFromLengths(codelengthcode, 7); 
        if(error) 
            return;
        size_t i = 0;
        size_t replength;
        while(i < HLIT + HDIST)
        {
            unsigned long code = huffmanDecodeSymbol(in, bp, codelengthcodetree, inlength); 
            if(error) 
                return;
            if(code <= 15)  { 
                if(i < HLIT) bitlen[i++] = code; 
                else bitlenD[i++ - HLIT] = code; 
            } //a length code
            else if(code == 16) //repeat previous
            {
                if(bp >> 3 >= inlength) //error, bit pointer jumps past memory
                { 
                    error = 50; 
                    return; 
                } 
                unsigned long newbits = readBitsFromStream(bp, in, 2);
                replength = (size_t) 3 + newbits;
                unsigned long value=0; //set value to the previous code
                if (i == 0) 
                    value = bitlen[0]; // BUG fixed here.
                else if((i - 1) < HLIT) 
                    value = bitlen[i - 1];
                else 
                    value = bitlenD[i - HLIT - 1];
                for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
                {
                    if(i >= HLIT + HDIST) //error: i is larger than the amount of codes
                    { 
                        error = 13; 
                        return; 
                    } 
                    if(i < HLIT) {
                        bitlen[i++] = value; 
                    }
                    else {
                        bitlenD[i++ - HLIT] = value;
                    }
                }
            }
            else if(code == 17) //repeat "0" 3-10 times
            {
                if(bp >> 3 >= inlength) //error, bit pointer jumps past memory
                { 
                    error = 50; 
                    return; 
                } 
                replength = 3 + readBitsFromStream(bp, in, 3);
                for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
                {
                    if (i >= HLIT + HDIST) //error: i is larger than the amount of codes
                    { 
                        error = 14; 
                        return; 
                    } 
                    if (i < HLIT) 
                        bitlen[i++] = 0; 
                    else 
                        bitlenD[i++ - HLIT] = 0;
                }
            }
            else if(code == 18) //repeat "0" 11-138 times
            {
                if(bp >> 3 >= inlength) //error, bit pointer jumps past memory
                { 
                    error = 50; 
                    return; 
                } 
                replength = 11 + readBitsFromStream(bp, in, 7);
                for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
                {
                    if(i >= HLIT + HDIST) //error: i is larger than the amount of codes
                    { 
                        error = 15; 
                        return; 
                    } 
                    if(i < HLIT) 
                        bitlen[i++] = 0; 
                    else 
                        bitlenD[i++ - HLIT] = 0;
                }
            }
            else //error: somehow an unexisting code appeared. This can never happen.
            { 
                error = 16; 
                return; 
            } 
        }
        if (bitlen[256] == 0) //the length of the end code 256 must be larger than 0
        { 
            error = 64; 
            return; 
        } 
        if (bitlen[0] > 288 ) // codes are too big and making crashiness // CQ
        { 
            error = 600; 
            return; 
        } 
        error = tree.makeFromLengths(bitlen, 15); 
        if(error) //now we've finally got HLIT and HDIST, so generate the code trees, and the function is done
            return; 
        error = treeD.makeFromLengths(bitlenD, 15); 
        if(error) 
            return;
    }

    int inflateHuffmanBlock(std::vector<unsigned char>& out, const unsigned char* in, size_t& bp, size_t& pos, size_t inlength, unsigned long btype) 
    {
        if(btype == 2) { 
            getTreeInflateDynamic(codetree, codetreeD, in, bp, inlength); 
            if(error)  {
                return 0; 
            } 
        } else {
            return 0; 
        }

        for(;;)
        {
            unsigned long code = huffmanDecodeSymbol(in, bp, codetree, inlength); 
            if(error) return 0;
            if(code == 256) 
                break;
            else if (code <= 255) //literal symbol
            {
                if(pos >= out.size()) 
                    out.resize((pos + 1) * 2); //reserve more room
                out[pos++] = (unsigned char)(code);
                //pos++;
            }
            else if (code >= 257 && code <= 285) //length code
            {
                size_t length = LENBASE[code - 257], numextrabits = LENEXTRA[code - 257];
                if ((bp >> 3) >= inlength) { 
                    error = 51; 
                    return 0; 
                } //error, bit pointer will jump past memory
                length += readBitsFromStream(bp, in, numextrabits);
                unsigned long codeD = huffmanDecodeSymbol(in, bp, codetreeD, inlength); 
                if (error) 
                    return 0;
                if (codeD > 29) { 
                    error = 18; 
                    return 0; 
                } //error: invalid dist code (30-31 are never used)
                unsigned long dist = DISTBASE[codeD], numextrabitsD = DISTEXTRA[codeD];
                if ((bp >> 3) >= inlength) { 
                    error = 51; 
                    return 0; 
                } //error, bit pointer will jump past memory
                dist += readBitsFromStream(bp, in, numextrabitsD);
                // CQ overflow in orig code fixed
                //uint64_t start = pos, back = start - dist; //backwards
                long start = pos;
                long back = pos ; 
                if (dist >= start) 
                    back=start; 
                else 
                    back=start-dist; //backwards
                if(pos + length >= out.size()) 
                    out.resize((pos + length) * 2); //reserve more room
                for(size_t i = 0; i < length; i++) { 
                    out[pos++]= out[back++];
                    if(back >= start)  {
                        if (dist >= start) 
                            back=start; 
                        else
                            back = start - dist; 
                    }
                }
            }
        } 
        // if we successfully decoded, deal with our trees
        // only worrying about tree #1 with list of codes
        int hufftablelist[256];
        for (int i=0; i < 256; i++)
            hufftablelist[i]=0;
        for (int i=0; i < codetree.tree2d.size()/2; i++){
            if (codetree.tree2d.at(i) < 256 ) {
                hufftablelist[codetree.tree2d.at(i)]=1;
            } 
        }
        int xmlct=0;
        int pngexect=0;
        int pngct=0;
        int elemct=0;
        for (int i=0; i < 256; i++) {
            if (i < 10) 
                xmlct+=hufftablelist[i]; 
            if (i > 108 && i < 120)
                pngexect+=hufftablelist[i]; 
            // second classifier for png
            if (i > 233 && i < 252) 
                pngct+=hufftablelist[i];
            elemct+=hufftablelist[i];
        }
        pngexect+=hufftablelist[255];
        int res=0;
        if (xmlct == 0) { 
            cout << "xml found";
            res=2;
        } else { 
            if (pngexect > 4 && pngexect < 8 && elemct < 110) {
                //cout << "png+ found";
                res=4;
            } else if (pngexect > 4 && pngct > 1) {
                //cout << "png+ found " << pngct ;
                res=4;
            } else if (pngexect > 4) {
                //cout << "exe found";
                res=3;
            } else {
                //cout << "png found";
                res=4;
            }
        }
        //cout << ", elements " << elemct;
        //cout << endl;
        error=0;
        return res;
    } // end inflateHuffmanBlock
}; // end of inflator

/**
decodeStuff: decodes a buffer in memory, into a raw data buffer. 
Slides a window to find possible DEFLATE blocks

in_image: pointer to the buffer of the source file in memory. To get it from a file on
disk, load it and store it in a memory buffer yourself first.
in_size: size of the input file in bytes.
return: 0 if success, not 0 if some error occured.
*/

int decodeStuff(std::vector<unsigned char>& out_image, const unsigned char* in_image, size_t in_size, const char *fragname)
{
    std::vector<unsigned char> idat; 
    idat.insert(idat.end(), &in_image[0], &in_image[in_size]);

    Inflator inflator;
    if (idat.size() < 2) { //error, size of data too small
        return 53;
    } 
    inflator.inflate(out_image, idat, 2);
    cout <<  "summary for " << fragname << ", " << inflator.blockcount << " deflate block(s) found"<< endl;
    if (inflator.blockcount > 0)  {
        if (inflator.xmlcount) {
            cout << "    xml " << inflator.xmlcount <<"/"<< inflator.blockcount<< " ";
            if (inflator.xmlcount==inflator.blockcount) 
                cout << "100%"<<endl;
            else
                cout << setprecision(2) << 100.0*(float) inflator.xmlcount / inflator.blockcount <<"%" << endl;
        }
        if (inflator.execount) {
            cout << "    exe " << inflator.execount <<"/"<< inflator.blockcount<< " ";
            if (inflator.execount==inflator.blockcount) 
                cout << "100%"<<endl;
            else 
                cout << setprecision(2) << 100.0*(float) inflator.execount / inflator.blockcount <<"%" << endl;
        }
        if (inflator.pngcount) {
            cout << "    png " << inflator.pngcount <<"/"<< inflator.blockcount<< " ";
            if (inflator.pngcount==inflator.blockcount) 
                cout << "100%"<<endl;
            else 
            cout << setprecision(2) << 100.0*(float) inflator.pngcount / inflator.blockcount <<"%" << endl;
        }
    }
    return inflator.error;
}

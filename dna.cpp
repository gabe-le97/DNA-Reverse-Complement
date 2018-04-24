/*
 * File: dna.cpp
 * Includes dna.h & main.cpp
 * Created by Gabe Le and Danny Nguyen
 *------------------------------------------------------------------------------
 * This file will take a fasta formatted file and give the reverse complement
 * while checking many rules. This is done through multiple functions &
 * constructors in the class DNA
 */
#include "dna.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
//------------------------------------------------------------------------------
/*
 * Default Constructor that sets the private variables equal to an empty string
 */
DNA::DNA(){
   hdr = "";
   seq = "";
}
//------------------------------------------------------------------------------
/* This constructor will give the private variables hdr & str the values of
 * sequence and header. It will check if the header and sequence contain valid
 * characters and if there are invalid characters it will throw a runtime error
 */
DNA::DNA(std::string header, std::string sequence) {
   hdr = header;
   seq = sequence;
// Searching the header for the > character to validate it
   if (header[0] != '>')  {
    throw std::runtime_error("The header is invalid!");
  }
  else {
// If true, hdr will hold the first line of the file (should be header)
    hdr = header;
  }
// Increments through the array to check if each character is valid for sequence
  for (int i = 0; i < (int) sequence.length(); i++){
    if (sequence[i] != 'A' && sequence[i] != 'T' && sequence[i] != 'C'
    && sequence[i] != 'G' && sequence[i] != 'N') {
    throw std::runtime_error("The Sequence contains an invalid character!");
    }
  }
// seq now holds the DNA string
   seq = sequence;
}
//------------------------------------------------------------------------------
/*
 * Getter method that returns the sequence to the function because it can't be
 * accessed directly (private)
 */
std::string DNA::getSequence(){
  return seq;
}
//------------------------------------------------------------------------------
/*
 * Setter method that returns the header to the function because it can't be
 * accessed directly (private)
 */
std::string DNA::getHeader(){
  return hdr;
}
//------------------------------------------------------------------------------
/*
 * Reverse method using a switch function that will return the value as a reverse
 * complement. The sequence will be traversed starting from the end to the
 * beginning
 */
DNA DNA::revcomp(){
  std::string revcomp = "";
// Decrement the seq length in order to traverse the string backwards
  for (int i = (int)seq.length() - 1; i >= 0; i--)
    switch(seq[i]){
      case 'A':
        revcomp += 'T'; break;
      case 'C':
        revcomp += 'G'; break;
      case 'G':
        revcomp += 'C'; break;
      case 'T':
        revcomp += 'A'; break;
      case 'N':
        revcomp += 'N'; break;
    }
// revcomp now holds the reverse complement of the sequence
// Runs hdr and revcomp(sequence) as arguments through the constructor
    DNA::DNA revcomp2 = DNA::DNA(hdr, revcomp);
    return revcomp2;
}
//------------------------------------------------------------------------------
/*
 * Function in order to find the starting point of a defined sequence & its reverse
 * complement. If query is not found, it will return std::string::npos (not found)
 */
size_t DNA::find(std::string query, size_t start){
// Looking for defined sequence in the sequence
  size_t nuc = seq.find(query, start);
// Running the constructor to check for any errors with the header and query
  DNA::DNA chromo = DNA::DNA(hdr, query);
  DNA::DNA mutate = chromo.revcomp();
// gene holds value of seq now
  std::string gene = mutate.getSequence();
// Looking for the query in gene now
  size_t electro = seq.find(gene, start);
  if (electro > nuc && nuc != std::string::npos){
// Returns the index of the query or its reverse complement
    return nuc;
  }
  else {
// Returns std::string::npos which is -1
    return electro;
  }
}
//------------------------------------------------------------------------------
/*
 * Constructor that calls for file input and checks if the header and sequence have valid characters
 */
DNA::DNA(std::ifstream &infile) {
// The file is pushed into a vector named riboflex
  std::vector <std::string> riboflex;
  std::string fileSequence = "";
  std::string fileHeader = "";
  std::string line;
// Now the file is being taken line by line and placed into line
    while (getline(infile, line)) {
      riboflex.push_back(line);
  }
//Increment through the string and make fileHeader a string of the first line of the file
  for (int i = 0; i < (int) riboflex.size(); i++){
    if (riboflex[i] == riboflex[0]){
      fileHeader += riboflex[i];
    }
// Now that riboflex is not 0, we will begin to fill the sequence string
    else
      fileSequence += riboflex[i];
    }
    if (fileHeader[0] != '>') {
      throw std::runtime_error("Invalid Header");
    }
// Increment through the sequence to check for invalid characters like before
      for (int i = 0; i < (int) fileSequence.length(); i++) {
        if (fileSequence[i] != 'A' && fileSequence[i] != 'T' &&
        fileSequence[i] != 'C' && fileSequence[i] != 'G' && fileSequence[i] != 'N') {
          throw std::runtime_error("Invalid sequence format.");
          }
        }
// private variables hdr & seq now hold a new value
  hdr = fileHeader;
  seq = fileSequence;
}
//------------------------------------------------------------------------------
/*
 * Prints out the fasta output as a string restricted to a column argument
 * By default column is 80
 */
std::string DNA::toFasta(int columns){
  std::string deox = hdr + '\n' + seq[0];
// Increment through the seq to make a newline where argument column specifies
  for (int i = 1; i < (int) seq.length(); i++){
    if (i % columns == 0) {
      deox += "\n";
    }
      deox += seq[i];
    }
  return deox += '\n';
 }
//------------------------------------------------------------------------------
/*
 * Function that tests for equality with d1 & d2 (DNA sequences) excluding headers
 */
bool operator == (DNA::DNA d1, DNA::DNA d2){
  std::string s = d1.getSequence();
  std::string t = d2.getSequence();
// Will return either true or false
  return (s == t);
}

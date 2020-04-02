/*
 * parser.cpp
 * 
 * Copyright 2019 Miquel Bernat Laporta i Granados <mlaportaigranados@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#include <iostream>
#include <string>
#include <fstream>
#include "include/exit_codes.h"
#include "include/parser.h"
using namespace std;

std::string s;
std::ifstream read;
std::ofstream write;

void open_files(ifstream& read,ofstream& write)//Definition for the open_files()
{
    string infile="prova.txt";//reading the file
    read.open(infile.c_str());//open the read file
    string outfile="resultats.txt";
    write.open(outfile.c_str());//open the write file
    process_data(read,write);//calling the process_data() function
    return;
}

void process_data(ifstream& read,ofstream& write)//Definition for the process_data()
{
	unsigned int number=0,position,flag1 = 0, flag2 = 0 , flag3 = 0;//I am using flags to check the different conditions we encounter
	if(read.is_open())
	{
	while(getline(read,s))//to read the whole file line by line
        {
		unsigned int i;
		string s1;
                position = s.size();
                flag2 = 0;
                flag3 = 0;
                for(i = 0 ;i<s.size(); i++)
                {
                        // this checks for the double quotes "
                        if(s[i] == '"' && flag1 == 0 && flag2 == 0 && number == 0)
                        {
  	                      flag3=1;
                              number = 1;
                              s1 += '"';
                              i++;
                        }
          		if(s[i] == '"' && flag1 == 0 && flag2 == 0 && number == 1)
              		{
                        	flag3 = 0;
                        	number = 0;
                        	s1 += '"';
                        	i++;
                        }

                        //this checks for Double-slash //
                        if(s[i] == '/' && s[i+1] == '/' && flag3 == 0 && flag1 == 0)
                        {
                                flag2 = 1;
                                position = i;
                        }

                        // This checks for slash-star /*
                        if(s[i] == '/' && s[i+1] == '*' && flag2 == 0 && flag3 == 0)
                        {
                                flag1 = 1;
                                i += 2;
                        }
                        if(s[i] == '*' && s[i+1] == '/' && flag2 == 0 && flag3 == 0)
                        {
                                flag1 = 0;
                                i += 2;
                                if(s[i] == '/' && s[i+1] == '*')
                                {
                                        flag1 = 1;
                                }
                        }
                        

                        if(flag1 == 0)
                        {
                                if(i < position)
                                {
                                        s1 += s[i];
                                }
                        }
                }
                write<<s1<<endl;
        
        }
}
else
{
error_message("Error in opening the file");
}
}

void close_files(ifstream& read,ofstream& write)//Definition for the close_file() function
{
read.close();//this is close the read file
write.close();//this is close the write file
}

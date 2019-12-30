/*
Parser program test, comment removal implementation.
*/

#include <iostream>
#include <string>
#include <fstream>

using namespace std;

string s;
ifstream read;
ofstream write;
//Prototype function declarations
void open_files(ifstream&,ofstream&);
void process_data(ifstream&,ofstream&);
void close_files(ifstream&,ofstream&);

void error(string msg)//this is for printing the error message.
{
	cerr<<msg<<endl;
	exit(EXIT_FAILURE);
}

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
error("Error in opening the file");
}
}

void close_files(ifstream& read,ofstream& write)//Definition for the close_file() function
{
read.close();//this is close the read file
write.close();//this is close the write file
}

int main(){
open_files(read,write);//calling the open_files() function
close_files(read,write);//calling the close_files() function
return 0;
}

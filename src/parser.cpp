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

#include "include/parser.h"
#include "include/structures.h"
#include "include/menu.h"
#include "include/integration.h"

//Generic non defined gravity constant for use with user desired systems
double gravity_constant;

/*1st pointer is for the file in which we check for comments,
and 2nd is the file in which we copy the code after removing comments  */
FILE *fp , *fp2;

void check_comment(char c)
{
    char d;

    if( c == '/')   // if the character starts with '/', it 'could' be a comment
    {
        if((d=fgetc(fp))=='*')   // if the next character we read is '*', it is the beginning of multiblock comment
            block_comment();  // pass control to function that handles multiblock comments

        else if( d == '/')   // else if the next character we read is '/', it is the beginning of single line comment
        {
            single_comment();// pass control to function that handles single line comment

        }
        else
        {
            // if both the cases fail, it is not a comment, so we add the character as it is in the new file.
            fputc(c,fp2);
            fputc(d,fp2);
        }
    }

        // again, if all above fails, we add the character as it is in the new file.
    else
        fputc(c,fp2);
}

void block_comment()
{

    char d,e;
    while((d=fgetc(fp))!=EOF)   // the block comment has started, read the character by character
    {
        /* keep reading the characters and do nothing,
        as they do not have to be copied into the new file (we are removing the comments)
        */
        if(d=='*')    // if the comment 'seems' like ending
        {
            e=fgetc(fp);  // check if it actually ends (block comments end with '*/')

            if(e=='/')  // if the comment 'has' ended, return from the function
                return;
        }
    }
}

void single_comment()
{
    char d;

    while((d=fgetc(fp))!=EOF)  // the single line comment has started, read the character by character
    {
        /* keep reading the characters and do nothing,
        as they do not have to be copied into the new file (we are removing the comments)
        */
        if(d=='\n')   // check if the comment ends (single comments end with '\n', or newline)
            return;  // if the comment 'has' ended, return from the function

    }
}

void parse_file(){
    char c;
    fp = fopen ("test.txt","r") ;   // open the first file in read mode
    fp2 = fopen ("mynewfile.txt","w") ;    // open the second file in write mode

    while((c=fgetc(fp))!=EOF)       // read the file character by character
        check_comment(c);   // check for each character if it seems like the beginning of a comment

    //  close both the files at the end of the program
    fclose(fp);
    fclose(fp2);
}

void parse_data(const std::string& filename){

    int num_bodies;
    std::ifstream fin;
    fin.open("input.txt");
    if (!fin) {
        error_message("Error in opening the required file");
    }

    std::vector<body> bodies;
    body temp;
    fin >> num_bodies >> gravity_constant;
    while(fin >> temp.name >> temp.mass >> temp.radius >> temp.location.x >> temp.location.y >> temp.location.z >>
    temp.velocity.x >> temp.velocity.y >> temp.velocity.z){
        bodies.push_back(temp);
    }

// now print the information you read in
    for (const auto& simulation_objects : bodies) {
        std::cout << simulation_objects.name << ' ' << simulation_objects.mass << ' ' << simulation_objects.radius << std::endl;
        std::cout << "Position vector:" << std::endl;
        std::cout << simulation_objects.location.x << ' ' << simulation_objects.location.y << ' ' << simulation_objects.location.z << std::endl;
        std::cout << "velocity vector:" << std::endl;
        std::cout << simulation_objects.velocity.x << ' ' << simulation_objects.velocity.y << ' ' << simulation_objects.velocity.z << std::endl;
    }
}

//C implementation of the parse file function: benchmark and refactor needed
/*
void initiateSystem(char* fileName){
    int i;
    FILE* fp1 = fopen(fileName,"r");
    std::vector<body> bodies;
    double timeStep;
    int num_bodies;
    fscanf(fp,"%lf%d%d",&gravity_constant,&num_bodies,&timeStep);

    masses = (double*)malloc(bodies*sizeof(double));
    positions = (vector*)malloc(bodies*sizeof(vector));
    velocities = (vector*)malloc(bodies*sizeof(vector));
    accelerations = (vector*)malloc(bodies*sizeof(vector));

    for(i=0; i < num_bodies; i++){
        fscanf(fp,"%lf",&masses[i]);
        fscanf(fp,"%lf%lf%lf",&positions[i].x,&positions[i].y,&positions[i].z);
        fscanf(fp,"%lf%lf%lf",&velocities[i].x,&velocities[i].y,&velocities[i].z);
    }

    fclose(fp);
}
*/
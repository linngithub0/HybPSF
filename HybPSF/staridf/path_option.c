/**
 * C program to check whether a directory exists or not.
 */

#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

//int isDirectoryExists(char *path);
int pathexist(char *path)
{ 
    int isDirectoryExists(char *path);
    //char path[100];
    int status;
    //printf("Enter directory path: ");
   //	scanf("%s", path);


    // Check if directory exists or not
    if (isDirectoryExists(path))
    {
        printf("Directory exists at path '%s'\n", path);
    }
    else
    {
        printf("Directory does not exists at path '%s', and making\n", path);
        status = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if(status==0){printf("mkdir success!\n");}
        else printf("failure\n");
    }

    return 0;
}



/**
 * Function to check whether a directory exists or not.
 * It returns 1 if given path is directory and  exists 
 * otherwise returns 0.
 */
int isDirectoryExists(char *path)
{
    struct stat stats;

    stat(path, &stats);

    // Check for file existence
    if (S_ISDIR(stats.st_mode))
        return 1;

    return 0;
}

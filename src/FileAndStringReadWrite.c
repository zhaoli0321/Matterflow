/*! \file  FileAndStringReadWrite.c
 *
 *  \brief Define some functions
 *
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2021--2022 by the CAEP-XTU team. All rights reserved.
 *---------------------------------------------------------------------------------
 */


#include "matterflow.h"
#include "matterflow_functs.h"


FILE * MyOpenFile(char * fileName, char * format)
{
	FILE * fp;
	fp = fopen(fileName, format);
	if (fp == NULL)
	{
		printf("cannot open file %s", fileName);
		ERR("\n");
	}
	return(fp);
}

bool SearchNextInFile(FILE * fp, char * searchStr)
{
	char chGet, *pChSrch;
	do
	{
		fscanf(fp, "%c", &chGet);
		if (chGet == (*searchStr))
		{
			pChSrch = searchStr;
			do
			{
				if ((*(++pChSrch)) == '\0')return(true);
				fscanf(fp, "%c", &chGet);
				if (chGet != (*pChSrch))break;
			} while (true);
			fseek(fp, (long)(searchStr - pChSrch), SEEK_CUR);
		}
	} while (!feof(fp));
	return(false);
}

bool ContainStr(char * strSource, char * searchStr)
{
	char * position, chGet, *pChSrch;
	position = strSource;
	do
	{
		chGet = *(position++);
		if (chGet == (*searchStr))
		{
			pChSrch = searchStr;
			do
			{
				if ((*(++pChSrch)) == '\0')return(true);
				chGet = *(position++);
				if (chGet != (*pChSrch))break;
			} while (true);
		}
	} while ((*position) != '\0');//注意这里修改了。
	return(false);
}

double StringToDouble(char * str)
{
	double value;
	sscanf(str, "%lf", &value);
	return(value);
}

int StringToInt(char * str)
{
	int value;
	sscanf(str, "%d", &value);
	return(value);
}



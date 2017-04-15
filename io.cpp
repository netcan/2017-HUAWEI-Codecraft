#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <sys/timeb.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>

#define MAX_LINE_LEN 55000

#define INLINE  static __inline
#ifdef _DEBUG
#define PRINT   printf
#else
#define PRINT(...)
#endif


INLINE void write_file(const bool cover, const char * const buff, const char * const filename);


void print_time(const char *head)
{
#ifdef _DEBUG
    struct timeb rawtime;
    struct tm * timeinfo;
    ftime(&rawtime);
    timeinfo = localtime(&rawtime.time);

    static int ms = rawtime.millitm;
    static unsigned long s = rawtime.time;
    int out_ms = rawtime.millitm - ms;
    unsigned long out_s = rawtime.time - s;
    ms = rawtime.millitm;
    s = rawtime.time;

    if (out_ms < 0)
    {
        out_ms += 1000;
        out_s -= 1;
    }
    printf("%s date/time is: %s \tused time is %lu s %d ms.\n", head, asctime(timeinfo), out_s, out_ms);
#endif
}

int read_file(char ** const buff, const unsigned int spec, const char * const filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
    	PRINT("Fail to open file %s, %s.\n", filename, strerror(errno));
        return 0;
    }
    PRINT("Open file %s OK.\n", filename);

    char line[MAX_LINE_LEN + 2];
    unsigned int cnt = 0;
    while ((cnt < spec) && !feof(fp))
    {
        line[0] = 0;
        if (fgets(line, MAX_LINE_LEN + 2, fp) == NULL)  continue;
        if (line[0] == 0)   continue;
        buff[cnt] = (char *)malloc(MAX_LINE_LEN + 2);
        strncpy(buff[cnt], line, MAX_LINE_LEN + 2 - 1);
        buff[cnt][MAX_LINE_LEN + 1] = 0;
        cnt++;
    }
    fclose(fp);
    PRINT("There are %d lines in file %s.\n", cnt, filename);

    return cnt;
}

void write_result(const char * const buff,const char * const filename)
{
	// 以覆盖的方式写入
    write_file(1, buff, filename);

}

void release_buff(char ** const buff, const int valid_item_num)
{
    for (int i = 0; i < valid_item_num; i++)
        free(buff[i]);
}

INLINE void write_file(const bool cover, const char * const buff, const char * const filename)
{
    if (buff == NULL)
        return;

    const char *write_type = cover ? "w" : "a";//1:覆盖写文件，0:追加写文件
    FILE *fp = fopen(filename, write_type);
    if (fp == NULL)
    {
        PRINT("Fail to open file %s, %s.\n", filename, strerror(errno));
        return;
    }
    PRINT("Open file %s OK.\n", filename);
    fputs(buff, fp);
    fputs("\n", fp);
    fclose(fp);
}


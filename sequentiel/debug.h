#ifndef DEBUG__H
#define DEBUG__H

#include <stdio.h>

/**
 * Debug header. Refer to
 * http://stackoverflow.com/questions/1644868/c-define-macro-for-debug-printing
 * to get more details
 */

/* Set it to 0 to disable debug macros */
#define DEBUG 1

#if __STDC_VERSION__ >= 199901L
	#define C99
#endif

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

/* C99 exclusive */
#ifdef C99
	#define PRINT(fmt, ...) \
			do { fprintf(stderr, fmt "\n", ##__VA_ARGS__); } while (0)
	#define LOG(fmt, ...) \
            //do { if (DEBUG) fprintf(stderr, "L: " fmt "\n", ##__VA_ARGS__); } while (0)
	#define ERROR(fmt, ...) \
            do { if (DEBUG) fprintf(stderr, ANSI_COLOR_RED "E: " fmt "\n" ANSI_COLOR_RESET, ##__VA_ARGS__); } while (0)
	#define SUCCESS(fmt, ...) \
			do { if (DEBUG) fprintf(stderr, ANSI_COLOR_GREEN "S: " fmt "\n" ANSI_COLOR_RESET, ##__VA_ARGS__); } while (0)
#endif


#endif

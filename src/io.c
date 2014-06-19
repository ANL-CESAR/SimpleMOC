#include "SimpleMOC_header.h"

// Prints program logo
void logo(int version)
{
	border_print();
	printf(
"              _____ _                 _      __  __  ____   _____ \n"
"             / ____(_)               | |    |  \\/  |/ __ \\ / ____|\n"
"            | (___  _ _ __ ___  _ __ | | ___| \\  / | |  | | |     \n"
"             \\___ \\| | '_ ` _ \\| '_ \\| |/ _ \\ |\\/| | |  | | |     \n"
"             ____) | | | | | | | |_) | |  __/ |  | | |__| | |____ \n"
"            |_____/|_|_| |_| |_| .__/|_|\\___|_|  |_|\\____/ \\_____|\n"
"                               | |                                \n"
"                               |_|                                \n"
	);
	border_print();
	printf("\n");
	center_print("Developed at", 79);
	center_print("The Massachusetts Institute of Technology", 79);
	center_print("and", 79);
	center_print("Argonne National Laboratory", 79);
	printf("\n");
	char v[100];
	sprintf(v, "Version: %d", version);
	center_print(v, 79);
	printf("\n");
	border_print();
}

// Prints Section titles in center of 80 char terminal
void center_print(const char *s, int width)
{
	int length = strlen(s);
	int i;
	for (i=0; i<=(width-length)/2; i++) {
		fputs(" ", stdout);
	}
	fputs(s, stdout);
	fputs("\n", stdout);
}

// Prints a border
void border_print(void)
{
	printf(
	"==================================================================="
	"=============\n");
}

// Prints comma separated integers - for ease of reading
void fancy_int( int a )
{
    if( a < 1000 )
        printf("%d\n",a);

    else if( a >= 1000 && a < 1000000 )
        printf("%d,%03d\n", a / 1000, a % 1000);

    else if( a >= 1000000 && a < 1000000000 )
        printf("%d,%03d,%03d\n", a / 1000000, (a % 1000000) / 1000, a % 1000 );

    else if( a >= 1000000000 )
        printf("%d,%03d,%03d,%03d\n",
               a / 1000000000,
               (a % 1000000000) / 1000000,
               (a % 1000000) / 1000,
               a % 1000 );
    else
        printf("%d\n",a);
}

// Prints out the summary of User input
void print_input_summary(Input input)
{
	printf("INPUT SUMMARY\n");
}

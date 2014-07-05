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
void print_input_summary(Input I)
{
	center_print("INPUT SUMMARY", 79);
	border_print();
	printf("%-35s%d\n", "x-axis assemblies:", I.x_assemblies);
	printf("%-35s%d\n", "y-axis assemblies:", I.y_assemblies);
	printf("%-35s%d\n", "coarse axial intervals:", I.cai);
	printf("%-35s%d\n", "fine axial intervals:", I.fai);
	printf("%-35s%d\n", "axial source expansion order:", I.axial_exp);
	printf("%-35s%.2lf\n", "radial ray separation:", I.radial_ray_sep);
	printf("%-35s%.2lf\n", "axial z-ray separation:", I.axial_z_sep);
	printf("%-35s%d\n", "azimuthal angles:", I.n_azimuthal);
	printf("%-35s%d\n", "polar angles:", I.n_polar_angles);
	printf("%-35s%d\n", "energy groups:", I.n_egroups);
	if(I.decompose == 0)
		printf("%-35s%s\n", "data decomposition:", "OFF");
	else
		printf("%-35s%s\n", "data decomposition:", "ON");
	printf("%-35s%d\n", "assemblies per axial sub-domain:", I.decomp_assemblies_ax);
	printf("%-35s%ld\n", "avg segments per track:", I.segments_per_track);
	printf("%-35s%.2lf\n", "assembly width:", I.assembly_width);
	printf("%-35s%.2lf\n", "reactor height:", I.height);
	printf("%-35s%ld\n", "Src regions per assembly:", I.num_source_regions_per_assembly);
	border_print();
}

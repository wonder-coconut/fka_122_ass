#include <stdio.h>
    
int
main()
{
    int b = 0;
    {
	b = 1;
    }
    printf("b: %i\n", b);
    return 0;
}

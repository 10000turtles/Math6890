#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    cout << "Hellow World :D" << endl;

    printf("Arguments: %d\n\n", argc);

    for (int i = 0; i < argc; i++)
    {
        printf("Entry %d: %s\n", i + 1, argv[i]);
    }

    return 0;
}
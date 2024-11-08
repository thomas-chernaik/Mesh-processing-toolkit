//
// Created by thomas on 08/11/24.
//

#include "manifoldtesting.h"

int main(int argc, char **argv)
{
    std::string filename;
    // if we have an argument, then we use that as the file name to convert, otherwise we prompt the user
    if (argc == 2)
    {
        filename = argv[1];
    } else if (argc == 1)
    {
        std::cout << "No filename was provided to convert, please enter one: " << std::endl;
        std::cin >> filename;
    } else
    {
        std::cerr << "Please enter at most one command line argument." << std::endl;
        return -1;
    }
    // check that the filename ends in .tri
    if (filename.size() <= 4)
    {
        std::cerr << "The filename is not for a .tri file" << std::endl;
        return -2;
    } else if (filename.substr(filename.size() - 4) != ".tri")
    {
        std::cerr << "The filename is not for a .tri file" << std::endl;
        return -2;
    }
    ManifoldTesting mt;
    mt.LoadMesh(filename);

    std::cout << "The mesh is manifold" << std::endl;

    return 0;
}
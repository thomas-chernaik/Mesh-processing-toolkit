//
// Created by thomas on 08/11/24.
//

#include <numeric>
#include "../src/ManifoldTester.h"
#include "../src/FaceIndex.h"


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
    // if the filename ends in .tri then we convert it to a .face file
    if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".tri")
    {
        std::cout << "Converting " << filename << " to a .face file" << std::endl;
        FaceIndex fi;
        fi.readFile(filename);
        // get the file name without the extension (last 4 characters)
        filename.resize(filename.size() - 4);
        filename = filename;
        fi.writeFile(filename);
        filename = filename + ".face";
    }

    // create a manifold tester and test the mesh
    ManifoldTester mt;
    mt.readFile(filename);
    mt.testManifold();

    // if the program hasn't exited then the mesh is manifold
    std::cout << "The mesh is manifold" << std::endl;

    // calculate the genus of the mesh
    auto genuses = mt.CalculateGenus();
    if (genuses.size() != 1)
    {
        for (int i = 0; i < genuses.size(); i++)
        {
            std::cout << "The genus of component " << i << " is " << genuses[i] << std::endl;
        }
    }
    std::cout << "The total genus of the mesh is " << std::accumulate(genuses.begin(), genuses.end(), 0) << std::endl;

    return 0;
}
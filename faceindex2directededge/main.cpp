#include "../src/directededge.h"

int main(int argc, char** argv)
{
    std::string filename;
    // if we have an argument, then we use that as the file name to convert, otherwise we prompt the user
    if (argc == 2)
    {
        filename = argv[1];
    }
    else if (argc == 1)
    {
        std::cout << "No filename was provided to convert, please enter one: " << std::endl;
        std::cin >> filename;
    }
    else
    {
        std::cerr << "Please enter at most one command line argument." << std::endl;
        return -1;
    }
    // check that the filename ends in .face
    if (filename.size() < 6)
    {
        std::cerr << "The filename is not for a .face file" << std::endl;
        return -2;
    }
    else if (filename.substr(filename.size() - 5) != ".face")
    {
        std::cerr << "The filename is not for a .face file" << std::endl;
        return -2;
    }
    std::cout << "Converting " << filename << " to a .diredge file" << std::endl;
    DirectedEdge f2de;
    f2de.readFile(filename);
    // get the file name without the extension (last 5 characters)
    filename.resize(filename.size() - 5);
    f2de.writeFile(filename);
    return 0;
}
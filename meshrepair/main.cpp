//
// Created by thomas on 08/11/24.
//

#include "../src/MeshRepair.h"
#include "../src/FaceIndex.h"

// enum for the file types
enum FileType
{
    TRI,
    FACE,
    DIREDGE
};


int main(int argc, char **argv)
{
    FileType fileType = FACE;
    std::string filenameWithoutExtension;
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

    // find the file type from the file name
    // get the location of the last . in the filename
    size_t lastDot = filename.find_last_of('.');
    // if there isn't a dot then the file type is unknown
    if (lastDot == std::string::npos)
    {
        std::cerr << "The filename does not have an extension" << std::endl;
        return -2;
    }
    // get the file extension
    std::string extension = filename.substr(lastDot);
    switch (extension.size())
    {
        case 4:
            if (extension == ".tri")
            {
                fileType = TRI;
            }
            break;
        case 5:
            if (extension == ".face")
            {
                fileType = FACE;
            }
            break;
        case 8:
            if (extension == ".diredge")
            {
                fileType = DIREDGE;
            }
            break;
        default:
            std::cerr << "The filename does not have a valid extension" << std::endl;
            return -2;
    }
    // get the filename without an extension
    filenameWithoutExtension = filename.substr(0, lastDot);


    // if the filename ends in .tri then we convert it to a .face file
    if (fileType == TRI)
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
    MeshRepair mr;
    // if the file is a .diredge file then we read it in as such
    if (fileType == DIREDGE)
    {
        mr.readFileDiredge(filename);
    } else
    {
        mr.readFile(filename);
    }
    // repair the mesh
    mr.repairMesh();
    // write the repaired mesh to the same filetype as the input
    filenameWithoutExtension = filenameWithoutExtension + "_repaired";
    switch (fileType)
    {
        case TRI:
            mr.writeRepairedMeshTri(filenameWithoutExtension);
            break;
        case FACE:
            mr.writeRepairedMeshFace(filenameWithoutExtension);
            break;
        case DIREDGE:
            mr.writeFile(filenameWithoutExtension);
            break;
    }

    return 0;
}
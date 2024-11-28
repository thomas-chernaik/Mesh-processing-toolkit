//
// Created by thomas on 28/11/24.
//

#include "progressBar.h"

int previousPercentage = 100;

void printProgress(int percentage)
{
    if (percentage == 0 && previousPercentage != 0)
    {
        std::cout << "Progress: ";
        // add on spaces for the remaining percentage / 5
        for (int i = percentage / 5; i < 19; i++)
        {
            std::cout << " ";
        }
        // add a pipe to show the end of the bar
        std::cout << "|";
    } else if (percentage == previousPercentage)
    {
        return;
    }
    else
    {
        // delete the spaces and pipe, and last progress character
        // delete the pipe
        std::cout << "\b";
        // delete the spaces
        for (int i = previousPercentage / 5; i < 19; i++)
        {
            std::cout << "\b";
        }
        // delete the last character
        std::cout << "\b";

        switch (percentage % 5)
        {
            case 1:
                std::cout << "▍";
                break;
            case 2:
                std::cout << "▌";
                break;
            case 3:
                std::cout << "▋";
                break;
            case 4:
                std::cout << "▊";
                break;
            case 0:
                std::cout << "▉ ";
                break;
        }
        if(percentage == 100)
        {
            std::cout << "\b";
        }
        std::cout.flush();
        // add on spaces for the remaining percentage / 5
        for (int i = percentage / 5; i < 19; i++)
        {
            std::cout << " ";
        }
        // add a pipe to show the end of the bar
        std::cout << "|";
        std::cout.flush();
    }
    if(percentage == 100)
    {
        std::cout << std::endl;
    }
    previousPercentage = percentage;
}

#include "xmlWriter.h"


int main()
{

    XMLAlignmentTree tree;
    tree.addAlignment("hello");
    tree.addAlignment("goodbye");
    tree.write("myFile.xml");
    return 0;
}

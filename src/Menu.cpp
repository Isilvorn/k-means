#include "../include/Menu.h"

Menu::Menu()
{
    lines = new list<string>;
}

Menu::~Menu()
{
    delete lines;
}

void Menu::addline(string strNewLine) {
    lines->push_back(strNewLine);
    return;
}

void Menu::load(string strFile) {
    ifstream fileInput;
    string   strLine;

    fileInput.open(strFile);
    if (fileInput.is_open()) {
        while (!fileInput.eof()) {
            getline(fileInput,strLine);
            addline(strLine);
        }
    }

    return;
}


void Menu::draw(unsigned int tab,unsigned int totalWidth) {
    list<string>::iterator litr;
    unsigned int           maxlength = totalWidth;       // variable for storing the length of the longest string in the menu
    int                    n = 0;                        // keeps track of the line number of the menu item

    // getting max line length

    litr=lines->begin();
     while (litr!=lines->end()) {
        if (litr->length() > maxlength) maxlength=litr->length();
        litr++;
    }

    for (unsigned int i=0; i<tab; i++)
        cout << "     ";
    for (unsigned int i=0; i<(maxlength+10); i++)
        cout << "*";
    cout << endl;

    litr=lines->begin();
    while (litr!=lines->end()) {
        for (unsigned int i=0; i<tab; i++)
            cout << "     ";
        cout << "* [" << left << setw(2) << n << "] " << *litr;
        for (unsigned int i=0; i<((maxlength+2)-litr->length()); i++)
            cout << " ";
        cout << "*" << endl;
        litr++;
        n++;
    }

    for (unsigned int i=0; i<tab; i++)
        cout << "     ";
    for (unsigned int i=0; i<(maxlength+10); i++)
        cout << "*";
    cout << endl;
    for (unsigned int i=0; i<tab; i++)
        cout << "     ";

    lastmax = maxlength;

    return;
}

int Menu::prompt(void) {
    list<string>::iterator litr;
    int                    n = 0;
    int                    nReturn; // return value

    litr=lines->begin();
     while (litr!=lines->end()) {
        n++;
        litr++;
    }

    cout << "> ";
    cin >> nReturn;

    if ((nReturn < 0) || (nReturn > (n-1))) {
        cout << "Invalid Entry!" << endl;
        return (-1);
    }

    return nReturn;
}

void Menu::addenda(string strDesc,bool boolBottomLine) {

    cout << "* " << strDesc;
    for (unsigned int i=0; i<((lastmax+7)-strDesc.length()); i++)
        cout << " ";
    cout << "*" << endl;

    if (boolBottomLine) {
        for (unsigned int i=0; i<(lastmax+10); i++)
            cout << "*";
        cout << endl;
    }

}

void Menu::addenda(string strDesc,double n,int decimals,int width,bool boolBottomLine) {

    cout << "* " << strDesc;
    cout << left << setprecision(decimals) << fixed << setw(width);
    cout << n;
    for (unsigned int i=0; i<((lastmax+7)-strDesc.length()-width); i++)
        cout << " ";
    cout << "*" << endl;

    if (boolBottomLine) {
        for (unsigned int i=0; i<(lastmax+10); i++)
            cout << "*";
        cout << endl;
    }

}

void Menu::addenda(string strDesc,int n,int width,bool boolBottomLine) {

    cout << "* " << strDesc;
    cout << left << setw(width);
    cout << n;
    for (unsigned int i=0; i<((lastmax+7)-strDesc.length()-width); i++)
        cout << " ";
    cout << "*" << endl;

    if (boolBottomLine) {
        for (unsigned int i=0; i<(lastmax+10); i++)
            cout << "*";
        cout << endl;
    }

}

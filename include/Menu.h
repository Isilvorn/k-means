#include <list>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#ifndef MENU_H
#define MENU_H


class Menu
{
    public:
        Menu();
        virtual ~Menu();

        void load(string);
        void addline(string);
        void draw(unsigned int=0,unsigned int=0);
        int  prompt(void);
        void addenda(string,bool=true);
        void addenda(string,double,int,int,bool=true);
        void addenda(string,int,int,bool=true);

    protected:
        list<string>           *lines;
    private:
        unsigned int           lastmax;
};

#endif // MENU_H

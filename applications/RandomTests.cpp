

#include <iostream>

class MyClass
{
public:
    int AddOne(int a)
    {
        return a + 1;
    }
    int SubOne(int a)
    {
        return a - 1;
    }

    using MCCallback = decltype(&MyClass::AddOne); // int (MyClass::*)(int);

    MCCallback GetCallback(bool DoAdd)
    {
        if (DoAdd)
            return &MyClass::AddOne;
        return &MyClass::SubOne;
    }
};


int main()
{

    MyClass myClass;
    MyClass::MCCallback Callback = myClass.GetCallback(false);

    std::cout << (myClass.*Callback)(1) << std::endl;
    return 0;
}

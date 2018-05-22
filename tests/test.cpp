struct A
{
    A(const int& foo)
        : bar(foo)
    {
    }

    int bar;
};

void myFunc(A someA)
{
}

int main()
{
    A a = A(int());
    myFunc(a);
}

/***********************************************************************
*
* see former changes at file stack.h.versioninfos.txt
*
***********************************************************************/
#ifndef _stack_h_
#define _stack_h_
#include <stddef.h>

template < class T > class Stack {
  protected:
        struct stackEl {        // element of stack
                T el_data;
                stackEl *next;  // pointer to next element
        };
        stackEl *first; // pointer to first element
        int nElements; // number of elements;

  public:
        Stack();
        int push(T);    // zero if no memory available
        int pop(T &);   // zero if empty stack
        int top(T &);
        int getNumberOfElements();
};

template < class T > Stack < T >::Stack()
{
        first = NULL;
        nElements = 0;
}

template < class T > int Stack < T >::push(T el)
{
        stackEl *temp = new stackEl;    // create new element

        if (!temp)
                return 0;       // no memory available
        temp->el_data = el;     // insert data
        temp->next = first;
        first = temp;   // temp is first element now
        nElements++;
        return 1;
}

template < class T > int Stack < T >::pop(T & el)
{
        if (!first)
                return 0;       // stack is empty
        el = first->el_data;    // get data
        stackEl *temp = first;

        first = temp->next;     // second element is first now
        delete temp;    // delete element

        nElements--;
        return 1;
}

template < class T > int Stack < T >::top(T & el)
{
        if (!first)
                return 0;       // stack is empty
        el = first->el_data;    // get data
        return 1;
}

template < class T > int Stack < T >::getNumberOfElements()
{
        return nElements;
}

template < class T > class extStack : public Stack < T >
// extended stack with some additional features
// e.g. access to all the elements
{
  public:
        int getElement(int elementNumber, T & el);
};

template < class T > int extStack < T >::getElement(int elementNumber, T & el)
{
        typename Stack<T>::stackEl * oneEl;
        if (elementNumber < 1)
                return 0;
        if (elementNumber > this->nElements)
                return 0;
        if (1 == elementNumber) {
                el = this->first->el_data;
                return 1;
        }
        oneEl = this->first;
        for (int n = 1; n < elementNumber; n++) {
                oneEl = oneEl->next;
        }
        el = oneEl->el_data;
        return 1;
};

#endif

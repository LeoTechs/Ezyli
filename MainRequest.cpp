//include the header file
#include "Request.h"
#include <iostream>
#include <string>
// main method
int main(void)
{
    //create a request
    Request request(1,"train","Paris","Lyon","France","France",1.0);
    //print the request
    request.print(request);
    return 0;
}

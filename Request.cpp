//include the header file
#include "Request.h"
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>

//constructor 
Request::Request(double id,std::string type, std::string depature, std::string arrival, std::string departure_city, std::string arrival_city, double departureTime)
{
    this->id = id;
    this->type = type;
    this->depature = depature;
    this->arrival = arrival;
    this->departure_city = departure_city;
    this->arrival_city = arrival_city;
    this->departureTime = departureTime;
}
//getters
double Request::get_id()
{
    return this->id;
}

std::string Request::get_type()
{
    return this->type;
}

std::string Request::get_depature()
{
    return this->depature;
}

std::string Request::get_arrival()
{
    return this->arrival;
}

std::string Request::get_departure_city()
{
    return this->departure_city;
}

std::string Request::get_arrival_city()
{
    return this->arrival_city;
}

double Request::get_departureTime()
{
    return this->departureTime;
}

double Request::get_arrivalTime()
{
    return this->arrivalTime;
}

//setters
void Request::set_id(double id)
{
    this->id = id;
}

void Request::set_type(std::string type)
{
    this->type = type;
}

void Request::set_depature(std::string depature)
{
    this->depature = depature;
}

void Request::set_arrival(std::string arrival)
{
    this->arrival = arrival;
}

void Request::set_departure_city(std::string departure_city)
{
    this->departure_city = departure_city;
}

void Request::set_arrival_city(std::string arrival_city)
{
    this->arrival_city = arrival_city;
}

void Request::set_departureTime(double departureTime)
{
    this->departureTime = departureTime;
}

void Request::set_arrivalTime(double arrivalTime)
{
    this->arrivalTime = arrivalTime;
}

//method to compare two requests
bool Request::compare_request_fn(Request request1, Request request2)
{
    if(request1.type == request2.type  && request1.departure_city == request2.departure_city && request1.arrival_city == request2.arrival_city && request1.departureTime == request2.departureTime)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//method to print a request
void Request::print_request_fn(Request request)
{
    std::cout << "-----------------------" << std::endl;
    std::cout << "REQUEST INFOS : " << std::endl;
     std::cout << "______________________-" << std::endl;
    std::cout << "Request id : " << request.id << std::endl;
    std::cout << "Request type : " << request.type << std::endl;
    std::cout << "Request depature : " << request.depature << std::endl;
    std::cout << "Request arrival : " << request.arrival << std::endl;
    std::cout << "Request departure_city : " << request.departure_city << std::endl;
    std::cout << "Request arrival_city : " << request.arrival_city << std::endl;
    std::cout << "Request departureTime : " << request.departureTime << std::endl;
    std::cout << "-----------------------" << std::endl;
}

//method to enter a request and return a request
Request Request::enter_request_fn()
{
    double id;
    std::string type;
    std::string depature;
    std::string arrival;
    std::string departure_city;
    std::string arrival_city;
    double departureTime;
    std::cout << "Enter the request id : " << std::endl;
    std::cin >> id;
    std::cout << "Enter the request type : " << std::endl;
    std::cin >> type;
    std::cout << "Enter the request depature : " << std::endl;
    std::cin >> depature;
    std::cout << "Enter the request arrival : " << std::endl;
    std::cin >> arrival;
    std::cout << "Enter the request departure_city : " << std::endl;
    std::cin >> departure_city;
    std::cout << "Enter the request arrival_city : " << std::endl;
    std::cin >> arrival_city;
    std::cout << "Enter the request departureTime : " << std::endl;
    std::cin >> departureTime;
    Request request(id,type,depature,arrival,departure_city,arrival_city,departureTime);
    return request;
}
//this function take a request accepted by driver and list of request and return a list of request that can travel with the driver
std::vector<Request> Request::list_request_accepted_fn(Request request, std::vector<Request> list_request)
{
    std::vector<Request> list_request_accepted;
    for(int i = 0; i < list_request.size(); i++)
    {
        if(request.compare_request_fn(request,list_request[i]))
        {
            list_request_accepted.push_back(list_request[i]);
        }
    }
    return list_request_accepted;
}

//the main method
int main()
{ 
    //initialisation of the list of request
    std::vector<Request> list_request;
    //create a request R1
    
    



}



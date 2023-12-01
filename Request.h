#include <string> // Add missing import for string
//definit la classe Request
#ifndef REQUEST_H
#define REQUEST_H

// class Request
class Request
{
public:
    double id;
    std::string type;
    std::string depature;
    std::string arrival;
    std::string departure_city; // Use fully qualified name for string
    std::string arrival_city; // Use fully qualified name for string
    double departureTime;
    double arrivalTime;
    // contructor 
    Request(double id,std::string type, std::string depature, std::string arrival, std::string departure_city, std::string arrival_city, double departureTime);

    //getters
    double get_id();
    std::string get_type();
    std::string get_depature();
    std::string get_arrival();
    std::string get_departure_city();
    std::string get_arrival_city();
    double get_departureTime();
    double get_arrivalTime();
    //setters
    void set_id(double id);
    void set_type(std::string type);
    void set_depature(std::string depature);
    void set_arrival(std::string arrival);
    void set_departure_city(std::string departure_city);
    void set_arrival_city(std::string arrival_city);
    void set_departureTime(double departureTime);
    void set_arrivalTime(double arrivalTime);

    //method to compare two requests
    bool compare_request_fn(Request request1, Request request2);

    //method to print a request
    void print_request_fn(Request request);

   //method to enter the request attributes and return a request
    Request enter_request_fn();
  //
  std::vector<Request> Request::list_request_accepted_fn(Request request, std::vector<Request> list_request);

};
#endif // REQUEST_H



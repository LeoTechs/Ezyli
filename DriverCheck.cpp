#include <iostream>
#include <map>
#include <string>
#include <vector>
// for vector sort
#include <algorithm>

// some constats

// indecies of event
#define EVENT_ID_INDEX 0
#define EVENT_REALITY_INDEX 1
#define EVENT_TYPE_INDEX 2
#define EVENT_PROPOSED_TIME_INDEX 3
#define EVENT_REAL_TIME_INDEX 4
#define EVENT_TRIP_TYPE_INDEX 5
#define EVENT_CITY_INDEX 6
#define EVENT_TRIP_ID_INDEX 7
// ride can have multi destinations
// then this index will be used to knows all ids of destinations
// separated by ','
#define EVENT_DEST_CITIES_INDEX 8

// trips indecies
#define TRIP_ID_INDEX 0
#define TRIP_PICKUP_INDEX 1
#define TRIP_DROPOFF_INDEX 2
#define TRIP_DEPARTURE_TIME_INDEX 3
#define WALK_TIME_TO_JOIN_TRIP_INDEX 4
#define WALK_TIME_AFTER_TRIP_INDEX 5

// for multiple pickups possibilities
#define PICKUP_POSSIBILITIES_POINT_ID 0
#define DROP_OFF_POSSIBILITIES_ID 0
#define WALK_TIME_TO_JOIN_PICKUP_INDEX 1
#define WALK_TIME_AFTER_DROP_OFF_INDEX 2

// some strings
#define EVENT_TYPE_PICKUP "p"
#define EVENT_TYPE_DROPOFF "d"
#define EVENT_TYPE_POSITION "pos"
#define TRIP_TYPE_RIDE "r"
#define TRIP_TYPE_SHARING "s"

// block cause
#define CURRENT_REQUESTS "CURRENT_REQUESTS"
#define OTHER_REQUESTS "OTHER_REQUESTS"

// define common functions
//  function that get time from a time matrix

// fill string
std::vector<std::string> split_string(std::string str, std::string delimiter)
{
    std::vector<std::string> result;
    size_t pos = 0;
    std::string token;
    while ((pos = str.find(delimiter)) != std::string::npos)
    {
        token = str.substr(0, pos);
        result.push_back(token);
        str.erase(0, pos + delimiter.length());
    }
    result.push_back(str);
    return result;
}

/**
 * Calculates the time between two locations based on a given time matrix.
 *
 * @param from_id The starting location id.
 * @param to_id The destination location id.
 * @param time_matrics The map containing time values between locations.
 * @return The time between the two locations.
 * @throws std::invalid_argument if the key is not found in the time matrix.
 */
double get_time_fn(std::string from_id, std::string to_id, std::map<std::string, double> time_matrics)
{
    // create key
    std::string key = from_id + '-' + to_id;
    // check if key exists in time_matrics
    // when the key is not found in the map, it returns the iterator to end()
    // we raise exception to avoid segmentation fault error
    if (time_matrics.find(key) == time_matrics.end())
    {
        throw std::invalid_argument("key " + key + " not found in time_matrics");
    }
    return time_matrics[from_id + '-' + to_id];
}

/**
 * Adjusts the distance trip time based on uncertainty coefficients.
 * It's mean apply the uncertainty coefficient to the trip time.
 * Because the driver can have some unexpected events that can increase the trip time.
 *
 * @param time The original trip time.
 * @param inter_travel_time_uncertainty_coef The uncertainty coefficient for intercity trips.
 * @param intra_travel_time_uncertainty_coef The uncertainty coefficient for intracity trips.
 * @param is_ic Flag indicating if the trip is intercity or intracity (default is false).
 * @return The adjusted trip time.
 */
double ajust_distanced_trip_time_fn(double time, double inter_travel_time_uncertainty_coef, double intra_travel_time_uncertainty_coef, bool is_ic = false)
{
    // if trip is intercity then use inter_travel_time_uncertainty_coef
    // else use intra_travel_time_uncertainty_coef
    if (is_ic)
    {
        return time * (1 + inter_travel_time_uncertainty_coef);
    }
    else
    {
        return time * (1 + intra_travel_time_uncertainty_coef);
    }
}

bool trip_is_ic(std::vector<std::string> ev1, std::vector<std::string> ev2)
{
    return ev1[EVENT_CITY_INDEX] != ev2[EVENT_CITY_INDEX];
}

// get the time to go from an event to another
double get_ev_time_fn(std::vector<std::string> from, std::vector<std::string> to, double inter_travel_time_uncertainty_coef, double intra_travel_time_uncertainty_coef, std::map<std::string, double> time_matrics)
{
    // print from id and to id
    std::cout << "from id index : " << from[EVENT_ID_INDEX] << std::endl;
    // if it's a pick and drop of a ride, then use the ride time
    if (from[EVENT_TYPE_INDEX] == EVENT_TYPE_PICKUP && to[EVENT_TYPE_INDEX] == EVENT_TYPE_DROPOFF && from[EVENT_TRIP_ID_INDEX] == to[EVENT_TRIP_ID_INDEX] && from[EVENT_TRIP_TYPE_INDEX] == TRIP_TYPE_RIDE)
    {
        // in this case, we take the value of the duration of the trip (the time requested by the client)
        double from_proposed_time = std::stod(from[EVENT_PROPOSED_TIME_INDEX]);
        double to_proposed_time = std::stod(to[EVENT_PROPOSED_TIME_INDEX]);
        double normal_time = to_proposed_time - from_proposed_time;
        normal_time = normal_time * (1 + intra_travel_time_uncertainty_coef);
        return normal_time;
    }

    // if it's a time from drop of a ride to another point
    // then we will take the time of the drop of drops that is too far to the next point
    if (from[EVENT_TYPE_INDEX] == EVENT_TYPE_DROPOFF && from[EVENT_TRIP_TYPE_INDEX] == TRIP_TYPE_RIDE)
    {

        std::string from_drops_cities = from[EVENT_DEST_CITIES_INDEX];
        // split the string by ',' to get all drops
        std::vector<std::string> drops_cities = split_string(from_drops_cities, ",");

        // if vector is empty then raise an error

        // if vector is empty then raise an error
        if (drops_cities.size() == 0)
        {
            throw std::invalid_argument("try to calculate time from drop of a ride to another point but the drop points are empty");
        }

        // worst possible from
        double worst_possible_time = 0;
        std::string worst_possible_id = drops_cities[0];
        for (size_t i = 0; i < drops_cities.size(); i++)
        {
            // get the time from drop to the next point
            double time = get_time_fn(drops_cities[i], to[EVENT_ID_INDEX], time_matrics);
            // if time is greater than worst_possible_time then update it
            if (time > worst_possible_time)
            {
                worst_possible_time = time;
                worst_possible_id = drops_cities[i];
            }
        }

        bool is_ic = worst_possible_id != to[EVENT_CITY_INDEX];

        // return the worst possible time
        return ajust_distanced_trip_time_fn(worst_possible_time, inter_travel_time_uncertainty_coef, intra_travel_time_uncertainty_coef, is_ic);
    }

    // if it's to a pick of an event, then use the max of normal time to reach and proposed time
    // if (to[EVENT_TYPE_INDEX] == EVENT_TYPE_PICKUP)
    // {
    //     double normal_time = get_time(from[EVENT_ID_INDEX], to[EVENT_ID_INDEX]);
    //     double proposed_time = std::stod(to[EVENT_PROPOSED_TIME_INDEX]);
    //     normal_time = ajust_distanced_trip_time(normal_time, from, to);
    //     return std::max(normal_time, proposed_time);
    // }

    // it's a normal trip. to = drop and from=pick or position
    double normal_time = get_time_fn(from[EVENT_ID_INDEX], to[EVENT_ID_INDEX], time_matrics);
    double proposed_time = std::stod(to[EVENT_PROPOSED_TIME_INDEX]);
    normal_time = ajust_distanced_trip_time_fn(normal_time, inter_travel_time_uncertainty_coef, intra_travel_time_uncertainty_coef, trip_is_ic(from, to));
    return normal_time;
}

void print_events_fn(std::string context, std::vector<std::vector<std::string>> events)
{
    // print all attributes of events
    std::cout << context << std::endl;
    // print each attribute like att_name : att_value
    for (size_t i = 0; i < events.size(); i++)
    {
        std::cout << "event_id : " << events[i][EVENT_ID_INDEX] << std::endl;
        std::cout << "event_reality : " << events[i][EVENT_REALITY_INDEX] << std::endl;
        std::cout << "event_type : " << events[i][EVENT_TYPE_INDEX] << std::endl;
        std::cout << "event_proposed_time : " << events[i][EVENT_PROPOSED_TIME_INDEX] << std::endl;
        std::cout << "event_real_time : " << events[i][EVENT_REAL_TIME_INDEX] << std::endl;
        std::cout << "event_trip_type : " << events[i][EVENT_TRIP_TYPE_INDEX] << std::endl;
        std::cout << "event_city : " << events[i][EVENT_CITY_INDEX] << std::endl;
        std::cout << "event_trip_id : " << events[i][EVENT_TRIP_ID_INDEX] << std::endl;

        // print a line to separate events
        std::cout << "-------------------------------------" << std::endl;
    }
}

/**
 * @brief This function fills the real time of the events.
 *
 * The real time is the time that the event will occur.
 *
 * @param end_index The index of the last event where we will stop. put -1 if you want to fill all the events.
 * @param events The list of events, where each event is represented by a vector of strings.
 * @param inter_city_travel_time_uncertainty_coef The coefficient for inter-city travel time uncertainty.
 * @param intra_city_travel_time_uncertainty_coef The coefficient for intra-city travel time uncertainty.
 * @param time_matrics The map of time metrics.
 * @param pick_duration The duration for pick events.
 * @param drop_duration The duration for drop events.
 * @param wait_time_between_trips The wait time between trips.
 * @return The updated list of events with the real time filled.
 *
 * @note Each pick event is followed by his drop event.
 * @note we fill it progressiveley
 * @note for each event we add the the travel time from the previous event
 * @note and if the previous was pick event we add the pick duration
 * @note and if the previous was drop event we add the drop duration and the wait time between trips
 * @note And finaly if calcalulated time is less then the proposed time we keep the proposed time
 *
 */
std::vector<std::vector<std::string>> fill_real_time_fn(int end_index, std::vector<std::vector<std::string>> events, double inter_city_travel_time_uncertainty_coef, double intra_city_travel_time_uncertainty_coef, std::map<std::string, double> time_matrics, double pick_duration, double drop_duration, double wait_time_between_trips)
{
    std::cout << "size " << std::endl;
    int size = events.size();

    // if the end index is -1 or greater than sier we take the size of the list
    int till = end_index == -1 ? size : end_index + 1;
    till = till > size ? size : till;

    // if events is empty return events
    if (size == 0)
    {
        return events;
    }

    // fill event [0]
    events[0][EVENT_REAL_TIME_INDEX] = events[0][EVENT_PROPOSED_TIME_INDEX];

    // loop on events and fill
    for (size_t i = 1; i < till; i++)
    {
        // get the time to travel from the previous event
        double time = get_ev_time_fn(events[i - 1], events[i], inter_city_travel_time_uncertainty_coef, intra_city_travel_time_uncertainty_coef, time_matrics);
        // add the time to the previous event time
        double calculated_time = std::stod(events[i - 1][EVENT_REAL_TIME_INDEX]) + time;
        // if the previous event was pick event add the pick duration
        if (events[i - 1][EVENT_TYPE_INDEX] == EVENT_TYPE_PICKUP)
        {
            calculated_time += pick_duration;
        }
        // if the previous event was drop event add the drop duration and the wait time between trips
        else if (events[i - 1][EVENT_TYPE_INDEX] == EVENT_TYPE_DROPOFF)
        {
            calculated_time += drop_duration + wait_time_between_trips;
        }

        double proposed_time = std::stod(events[i][EVENT_PROPOSED_TIME_INDEX]);
        // if the calculated time is less then the proposed time we keep the proposed time
        std::cout << "calculated_time " << calculated_time << "proposed_time" << proposed_time << std::endl;
        if (calculated_time < proposed_time)
        {
            events[i][EVENT_REAL_TIME_INDEX] = std::to_string(proposed_time);
        }
        // else we keep the calculated time
        else
        {
            events[i][EVENT_REAL_TIME_INDEX] = std::to_string(calculated_time);
        }
    }
    return events;
}

/**
 * @brief Represents a vector of vectors of strings.
 *
 * This data structure is used to store a collection of events, where each event is represented as a vector of strings.
 * The outer vector represents the collection of events, while the inner vectors represent individual events.
 * Each inner vector contains a sequence of strings that describe the properties of the corresponding event.
 *
 * @note The elements in the inner vectors can be accessed using the subscript operator [].
 * @note The elements in the outer vector can be accessed using the subscript operator [].
 *
 * @see std::vector
 * @see std::string
 *
 * @note For the Logic of the function
 * @note 1. The position event must be the first event in the list
 * @note 2. The pick event must follow the drop event of the same trip
 * @note 3. The normal events pick events are ordered by the EVENT_PROPOSED_TIME_INDEX index (the time when the event is supposed to occur)
 */
std::vector<std::vector<std::string>> order_events_fn(std::vector<std::vector<std::string>> events)
{

    std::vector<std::vector<std::string>> ordered_events;

    // we find the index of the position event and add it to as the first item of the list
    auto it = std::find_if(events.begin(), events.end(), [](const std::vector<std::string> &event)
                           { return event[EVENT_TYPE_INDEX] == EVENT_TYPE_POSITION; });

    // when item is found, it will not be equal to end()
    //  end() point to the last element + 1
    if (it != events.end())
    {
        int position_index = it - events.begin();
        ordered_events.push_back(events[position_index]);
    }

    // we find all picks and order them according to the EVENT_PROPOSED_TIME_INDEX index
    std::vector<std::vector<std::string>> picks;
    for (auto event : events)
    {
        if (event[EVENT_TYPE_INDEX] == EVENT_TYPE_PICKUP)
        {
            picks.push_back(event);
        }
    }

    std::sort(picks.begin(), picks.end(), [](const std::vector<std::string> &event1, const std::vector<std::string> &event2)
              { return event1[EVENT_PROPOSED_TIME_INDEX] < event2[EVENT_PROPOSED_TIME_INDEX]; });

    // we loop in the list and insert drop element after the pick element
    // to find where is the drop of the pick, we use the EVENT_TRIP_ID_INDEX (it must be the same)
    for (auto pick : picks)
    {
        ordered_events.push_back(pick);
        for (auto event : events)
        {
            if (event[EVENT_TYPE_INDEX] == EVENT_TYPE_DROPOFF && event[EVENT_TRIP_ID_INDEX] == pick[EVENT_TRIP_ID_INDEX])
            {
                ordered_events.push_back(event);
            }
        }
    }

    return ordered_events;
}

class DriverOverlapsChecker
{

public:
    // map for time matrics
    std::map<std::string, double> time_matrics;
    // list of trips
    std::vector<std::vector<std::string>> TRIPS;
    // list of requests that driver can accept after
    std::vector<std::map<std::string, std::vector<std::string>>> L;
    std::map<std::string, std::vector<std::string>> R;
    double inter_city_travel_time_uncertainty_coef;
    double intra_city_travel_time_uncertainty_coef;
    double ride_time_uncertainty_coef;
    double pick_duration;
    double drop_duration;
    double wait_time_between_trips;
    double max_late_time;

    // constructor
    DriverOverlapsChecker()
    {
    }
    void initTravel(std::map<std::string, std::vector<std::string>> R, std::vector<std::map<std::string, std::vector<std::string>>> L)
    {
        this->R = R;
        this->L = L;
    }

    void init(std::map<std::string, std::vector<std::string>> R, std::map<std::string, double> time_matrics, std::vector<std::map<std::string, std::vector<std::string>>> L, std::vector<std::vector<std::string>> TRIPS, double inter_city_travel_time_uncertainty_coef, double intra_city_travel_time_uncertainty_coef, double ride_time_uncertainty_coef, double pick_duration, double drop_duration, double wait_time_between_trips, double max_late_time)
    {
        this->time_matrics = time_matrics;
        this->TRIPS = TRIPS;
        this->L = L;
        this->R = R;
        this->inter_city_travel_time_uncertainty_coef = inter_city_travel_time_uncertainty_coef;
        this->intra_city_travel_time_uncertainty_coef = intra_city_travel_time_uncertainty_coef;
        this->ride_time_uncertainty_coef = ride_time_uncertainty_coef;
        this->pick_duration = pick_duration;
        this->drop_duration = drop_duration;
        this->wait_time_between_trips = wait_time_between_trips;
        this->max_late_time = max_late_time;
    }

    int event_index(std::vector<std::string> event, std::vector<std::vector<std::string>> events)
    {
        auto it = std::find_if(events.begin(), events.end(), [event](const std::vector<std::string> &event2)
                               { return event2[EVENT_ID_INDEX] == event[EVENT_ID_INDEX]; });

        if (it == events.end())
        {
            // throw exception
            throw std::invalid_argument("event not found");
        }

        return std::distance(events.begin(), it);
    }

    std::vector<std::string> exclude_event_ids(std::vector<std::map<std::string, std::vector<std::string>>> LR)
    {
        // build new list called events_with_r from TRIPS
        // print started
        std::cout << "started" << std::endl;

        std::vector<std::string> result;
        std::vector<std::vector<std::string>> events_with_r(TRIPS);
        std::cout << "before pick" << std::endl;
        std::vector<std::string> pick_r = R["pick"];
        std::vector<std::string> drop_r = R["drop"];
        std::cout << "after pick" << std::endl;

        events_with_r.push_back(pick_r);
        events_with_r.push_back(drop_r);

        std::cout << "after PUSH BACK" << std::endl;
        print_events_fn("events_with_r before order", events_with_r);

        events_with_r = order_events_fn(events_with_r);

        print_events_fn("events_with_r afer order", events_with_r);

        std::cout << "after ORDER" << std::endl;

        int pick_r_index = event_index(pick_r, events_with_r);

        std::cout << "after PICK R INDEX" << std::endl;

        // fill
        events_with_r = fill_real_time_fn(pick_r_index, events_with_r, this->inter_city_travel_time_uncertainty_coef, this->intra_city_travel_time_uncertainty_coef, this->time_matrics, this->pick_duration, this->drop_duration, this->wait_time_between_trips);
        std::cout << "pick_r_index" << pick_r_index << std::endl;
        double r_real_departure_without_r;

        std::cout << "pick_ri_index2" << events_with_r[pick_r_index][4] << std::endl;
        double r_real_departure = std::stod(events_with_r[pick_r_index][EVENT_REAL_TIME_INDEX]);

        double proposed_duration_r = std::stod(pick_r[EVENT_PROPOSED_TIME_INDEX]);

        int L_size = LR.size();
        for (int i = 0; i < L_size; i++)
        {

            std::map<std::string, std::vector<std::string>> Ri = LR[i];
            std::vector<std::vector<std::string>> events_with_r_ri(events_with_r);

            std::vector<std::string> pick_ri = Ri["pick"];
            std::vector<std::string> drop_ri = Ri["drop"];

            events_with_r_ri.push_back(pick_ri);
            events_with_r_ri.push_back(drop_ri);
            events_with_r_ri = order_events_fn(events_with_r_ri);

            pick_r_index = event_index(pick_r, events_with_r_ri);

            int pick_ri_index = event_index(pick_ri, events_with_r_ri);

            events_with_r_ri = fill_real_time_fn(events_with_r_ri.size(), events_with_r_ri, this->inter_city_travel_time_uncertainty_coef, this->intra_city_travel_time_uncertainty_coef, this->time_matrics, this->pick_duration, this->drop_duration, this->wait_time_between_trips);
            print_events_fn("events_with_r_ri", events_with_r_ri);
            std::cout << "[pick_ri_index][EVENT_REAL_TIME_INDEX]" << events_with_r_ri[pick_ri_index][EVENT_REAL_TIME_INDEX] << std::endl;
            std::cout << "pick_ri_index" << pick_ri_index << std::endl;

            double ri_real_departure_without_r = std::stod(events_with_r_ri[pick_ri_index][EVENT_REAL_TIME_INDEX]);
            double ri_real_departure_impacted_by_r = std::stod(events_with_r_ri[pick_ri_index][EVENT_REAL_TIME_INDEX]);
            double r_real_departure_impacted_by_ri = std::stod(events_with_r_ri[pick_r_index][EVENT_REAL_TIME_INDEX]);

            double ri_real_departure = std::stod(events_with_r_ri[pick_ri_index][EVENT_REAL_TIME_INDEX]);
            double proposed_duration_ri = std::stod(pick_ri[EVENT_PROPOSED_TIME_INDEX]);

            std::cout << "proposed_duration_r " << proposed_duration_r << " proposed_duration_ri " << proposed_duration_ri << std::endl;
            std::cout << "r_real_departure_impacted_by_ri : " << r_real_departure_impacted_by_ri << " r_real_departure : " << r_real_departure << std::endl;

            if (proposed_duration_r <= proposed_duration_ri)
            {
                std::vector<std::vector<std::string>> time_trip_without_r(TRIPS);
                time_trip_without_r.push_back(pick_ri);
                time_trip_without_r.push_back(drop_ri);
                time_trip_without_r = order_events_fn(time_trip_without_r);
                int ri_pick_index = event_index(pick_ri, time_trip_without_r);
                time_trip_without_r = fill_real_time_fn(ri_pick_index, time_trip_without_r, this->inter_city_travel_time_uncertainty_coef, this->intra_city_travel_time_uncertainty_coef, this->time_matrics, this->pick_duration, this->drop_duration, this->wait_time_between_trips);
                double ri_real_departure_without_r = std::stod(time_trip_without_r[ri_pick_index][EVENT_REAL_TIME_INDEX]);

                std::cout << "ri_real_departure_without_r " << ri_real_departure_without_r << " ri_real_departure_impacted_by_r " << ri_real_departure_impacted_by_r << std::endl;
                if (ri_real_departure_without_r < ri_real_departure_impacted_by_r)
                {
                    double diff = ri_real_departure_impacted_by_r - ri_real_departure_without_r;
                    std::cout << "diff : " << diff << std::endl;
                    if (diff > max_late_time)
                    {
                        // print diff
                        std::cout << "diff : " << diff << std::endl;
                        result.push_back(pick_ri[EVENT_TRIP_ID_INDEX]);
                    }
                }
            }
            else
            {
                std::cout << "r_real_departure_impacted_by_ri < r_real_departure " << r_real_departure << r_real_departure_impacted_by_ri << std::endl;
                if (r_real_departure_impacted_by_ri < r_real_departure)
                {
                    // log("r_real_departure_impacted_by_ri < r_real_departure")
                    double diff = r_real_departure - r_real_departure_impacted_by_ri;
                    if (diff > max_late_time)
                    {
                        std::cout << "diff : " << diff << std::endl;
                        result.push_back(pick_ri[EVENT_TRIP_ID_INDEX]);
                    }
                }
            }
        }

        return result;
    }
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>> run()
    {

        std::vector<std::string> can_travel_with_r;
        std::vector<std::string> exclude;
        std::vector<std::string> no_exclude;

        bool R_is_sharing = R["pick"][EVENT_TRIP_TYPE_INDEX] == TRIP_TYPE_SHARING;

        if (R_is_sharing)
        {
            // call canTravelWithR_fn with L
            std::vector<std::map<std::string, std::vector<std::string>>> result = canTravelWithR_fn(R, L);

            // extract trips ids from result
            for (size_t i = 0; i < result.size(); i++)
            {
                std::string trip_id = result[i]["pick"][EVENT_TRIP_ID_INDEX];
                can_travel_with_r.push_back(trip_id);

                // remove elements from L that have the same trip id
                L.erase(std::remove_if(L.begin(), L.end(), [trip_id](std::map<std::string, std::vector<std::string>> &event)
                                       { return event["pick"][EVENT_TRIP_ID_INDEX] == trip_id; }),
                        L.end());
            }
            // save in can_travel_with_r
            // remove trips ids from L

            // Loop in L and remove all elements taht is with the same trip id

            // for (size_t i = 0; i < can_travel_with_r.size(); i++)
            // {
            //     // L.erase(std::remove(L.begin(), L.end(), result[i]), L.end());

            //     // remove all elements from L with the same trip id
            //     L.erase(std::remove_if(L.begin(), L.end(), [can_travel_with_r, i](const std::map<std::string, std::vector<std::string>> &event)
            //                            { return event["pick"][EVENT_TRIP_ID_INDEX] == can_travel_with_r[i]; }),
            //             L.end());
            // }
        }

        // loop in L check if L[i] is sharing else add to no_exclude
        for (size_t i = 0; i < L.size(); i++)
        {
            if (L[i]["pick"][EVENT_TRIP_TYPE_INDEX] == TRIP_TYPE_SHARING)
            {
                exclude.push_back(L[i]["pick"][EVENT_TRIP_ID_INDEX]);
            }
            else
            {
                no_exclude.push_back(L[i]["pick"][EVENT_TRIP_ID_INDEX]);
            }
        }
        // return can_travel_with_r, exclude, no_exclude
        return std::make_tuple(can_travel_with_r, exclude, no_exclude);
    }

    // this function take a request accepted by driver and list of request and return a list of request that can travel with the driver
    //  the function take a request accepted by driver and list of request and return a list of request that can travel with the driver
    std::vector<std::map<std::string, std::vector<std::string>>> canTravelWithR_fn(std::map<std::string, std::vector<std::string>> R, std::vector<std::map<std::string, std::vector<std::string>>> L)
    {
        // init result
        std::vector<std::map<std::string, std::vector<std::string>>> result;
        // loop on L
        for (size_t i = 0; i < L.size(); i++)
        {
            // get current request
            std::map<std::string, std::vector<std::string>> current_request = L[i];
            // check if current request can travel with R
            if (travelCompare_fn(R, current_request))
            {
                // add current request to result
                result.push_back(current_request);
            }
        }
        return result;
    };

    // this function take a request accepted by driver and a request from list of request and return true if the driver can travel with the request
    bool travelCompare_fn(std::map<std::string, std::vector<std::string>> R, std::map<std::string, std::vector<std::string>> C)
    {
        // if the request client is a pickup and the request accepted by the driver is a pickup and the two requests have the same departure time, the two have the same departure and arrival place
        bool r_is_sharing = R["pick"][EVENT_TRIP_TYPE_INDEX] == TRIP_TYPE_SHARING;
        bool c_is_sharing = C["pick"][EVENT_TRIP_TYPE_INDEX] == TRIP_TYPE_SHARING;
        bool c_and_r_is_sharing = c_is_sharing && r_is_sharing;
        bool c_and_r_has_same_departure_time = C["pick"][EVENT_PROPOSED_TIME_INDEX] == R["pick"][EVENT_PROPOSED_TIME_INDEX];
        bool c_and_r_has_same_departure_city = C["pick"][EVENT_CITY_INDEX] == R["pick"][EVENT_CITY_INDEX];
        bool c_and_r_has_same_arrival_city = C["drop"][EVENT_CITY_INDEX] == R["drop"][EVENT_CITY_INDEX];
        //print R is sharing
        std::cout << "R is sharing " << R["pick"][EVENT_TRIP_TYPE_INDEX] << std::endl;
        std::cout << "c_is_sharing " << c_is_sharing << std::endl;
        std::cout << "c_and_r_has_same_departure_time " << c_and_r_has_same_departure_time << std::endl;
        std::cout << "c_and_r_has_same_departure_city " << c_and_r_has_same_departure_city << std::endl;
        std::cout << "c_and_r_has_same_arrival_city " << c_and_r_has_same_arrival_city << std::endl;
        // print C arrival_city
        std::cout << "C arrival_city" << C["drop"][EVENT_CITY_INDEX] << std::endl;
        if (C["pick"][EVENT_TYPE_INDEX] == EVENT_TYPE_PICKUP && R["pick"][EVENT_TYPE_INDEX] == EVENT_TYPE_PICKUP && c_and_r_has_same_departure_time && c_and_r_has_same_departure_city && c_and_r_has_same_arrival_city&&c_and_r_is_sharing)
        {
            return true;
        }

        return false;
    }
};

int main(void)
{
    // Create an instance of DriverOverlapsChecker
    DriverOverlapsChecker checker;
    std::map<std::string, std::vector<std::string>> R;
    std::map<std::string, double> time_matrics;
    std::vector<std::map<std::string, std::vector<std::string>>> L;
    std::vector<std::vector<std::string>> TRIPS;
    double inter_city_travel_time_uncertainty_coef = 0.5;
    double intra_city_travel_time_uncertainty_coef = 0.5;
    double ride_time_uncertainty_coef = 0.1;
    double pick_duration = 5 * 60;
    double drop_duration = 5 * 60;
    double wait_time_between_trips = 5 * 60;
    double max_late_time = 1.0;
    // init timematrix
    time_matrics["poste-biyemassi"] = 600;
    time_matrics["biyemassi-nkolbisson"] = 20 * 60;
    time_matrics["nkolbisson-etougebe"] = 15 * 60;
    time_matrics["etougebe-poste"] = 15 * 20;
    time_matrics["poste-kondengui"] = 15 * 60;
    time_matrics["kondengui-mendong"] = 15 * 60;
    time_matrics["mendong-poste"] = 15 * 60;
    time_matrics["biyemassi-akwa"] = 4 * 60 * 60;
    time_matrics["akwa-kamkop"] = 4 * 60 * 60;
    time_matrics["kamkop-bonaberi"] = 5 * 60 * 60;
    time_matrics["bonaberi-biyemassi"] = 4 * 60 * 60;
    time_matrics["biyemassi-bonaberi"] = 4 * 60 * 60;
    time_matrics["bonaberi-baleveng"] = 5 * 60 * 60;
    time_matrics["baleveng-bonaberi"] = 5 * 60 * 60;
    time_matrics["bonaberi-olembe"] = 4 * 60 * 60;
    time_matrics["olembe-bonaberi"] = 4 * 60 * 60;
    time_matrics["bonaberi-bonendale"] =  20 * 60 ;
    time_matrics["bonendale-bonaberi"] =  20 * 60 ;
    time_matrics["bonaberi-yassa"] =      45* 60 ;
    time_matrics["yassa-bonaberi"] =      45* 60 ;
    time_matrics["bonaberi-pk14"] =       45* 60 ;
    time_matrics["pk14-bonaberi"] =       45* 60 ;
    //define time matrics poste-bepanda and bepanda-akwa
    time_matrics["poste-bepanda"] = 4 * 60 * 60;
    time_matrics["bepanda-akwa"] =      15 * 60;
    time_matrics["bepanda-olembe"] = 4* 60 * 60;
    time_matrics["olembe-yassa"] =  4 * 60 * 60;
    //define time matrics bepanda-bonaberi
    time_matrics["bepanda-bonaberi"] =  20 * 60;
    //define time matrics bepanda-bonendale
    time_matrics["bepanda-bonendale"] =  20 * 60;
    //define time matrics bonendale-pk14
    time_matrics["bonendale-pk14"] =  45 * 60;


    R["pick"] = {
        "poste",
        "1",
        EVENT_TYPE_PICKUP,
        "0",
        "0",
        "s",
        "yaounde",
        "r",
        ""};

    R["drop"] = {
        "bepanda",
        "1",
        EVENT_TYPE_DROPOFF,
        "10000",
        "0",
        "s",
        "douala",
        "r",
        ""};

    std::vector<std::string> pick_1 = {
        "nkolbisson",
        "1",
        EVENT_TYPE_PICKUP,
        "0",
        "0",
        "s",
        "yaounde",
        "r1",
        ""};

    // drop1

    std::vector<std::string> drop_1 = {
        "etougebe",
        "1",
        EVENT_TYPE_DROPOFF,
        "10000",
        "0",
        "s",
        "douala",
        "r1",
        ""};

    // pick 2
    std::vector<std::string> pick_2 = {
        "akwa",
        "1",
        EVENT_TYPE_PICKUP,
        "0",
        "0",
        "s",
        "douala",
        "r2",
        ""};

    // drop2
    std::vector<std::string> drop_2 = {
        "kamkop",
        "1",
        EVENT_TYPE_DROPOFF,
        "4000",
        "0",
        "d",
        "bafoussam",
        "r2",
        ""};

    //pick 3
    std::vector<std::string> pick_3 = {
        "bonaberi",
        "1",
        EVENT_TYPE_PICKUP,
        "0",
        "0",
        "s",
        "douala",
        "r3",
        ""};
    
    // drop3
    std::vector<std::string> drop_3 = {
        "baleveng",
        "1",
        EVENT_TYPE_DROPOFF,
        "4000",
        "0",
        "d",
        "dschang",
        "r3",
        ""};

    //pick 4
    std::vector<std::string> pick_4 = {
        "olembe",
        "1",
        EVENT_TYPE_PICKUP,
        "0",
        "0",
        "d",
        "yaounde",
        "r4",
        ""};

    // drop4
    std::vector<std::string> drop_4 = {
        "yassa",
        "1",
        EVENT_TYPE_DROPOFF,
        "4000",
        "0",
        "d",
        "douala",
        "r4",
        ""};

    //pick 5
    std::vector<std::string> pick_5 = {
        "bonendale",
        "1",
        EVENT_TYPE_PICKUP,
        "0",
        "0",
        "d",
        "douala",
        "r5",
        ""};
    // drop5
    std::vector<std::string> drop_5 = {
        "pk14",
        "1",
        EVENT_TYPE_DROPOFF,
        "4000",
        "0",
        "d",
        "douala",
        "r5",
        ""};
    
    // add to L
    std::map<std::string, std::vector<std::string>> R1;
    std::map<std::string, std::vector<std::string>> R2;
    std::map<std::string, std::vector<std::string>> R3;
    std::map<std::string, std::vector<std::string>> R4;
    std::map<std::string, std::vector<std::string>> R5;


    R1["pick"] = pick_1;
    R1["drop"] = drop_1;
    R2["pick"] = pick_2;
    R2["drop"] = drop_2;
    R3["pick"] = pick_3;
    R3["drop"] = drop_3;
    R4["pick"] = pick_4;
    R4["drop"] = drop_4;
    R5["pick"] = pick_5;
    R5["drop"] = drop_5;



    // add R1 to L
    L.push_back(R1);
    // add R2 to L
    L.push_back(R2);
    //push R3 to L
    L.push_back(R3);
    //push R4 to L
    L.push_back(R4);
    //push R5 to L
    L.push_back(R5);


    // call the init method
    checker.init(R, time_matrics, L, TRIPS, inter_city_travel_time_uncertainty_coef, intra_city_travel_time_uncertainty_coef, ride_time_uncertainty_coef, pick_duration, drop_duration, wait_time_between_trips, max_late_time);
    // call run
    std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::string>> result = checker.run();
    // get result
    std::vector<std::string> can_travel_with_r = std::get<0>(result);
    std::vector<std::string> exclude = std::get<1>(result);
    std::vector<std::string> no_exclude = std::get<2>(result);

   
    // print result
    std::cout << "can_travel_with_r" << std::endl;
    for (size_t i = 0; i < can_travel_with_r.size(); i++)
    {
        std::cout << can_travel_with_r[i] << std::endl;
    }
    std::cout << "exclude" << std::endl;
    for (size_t i = 0; i < exclude.size(); i++)
    {
        std::cout << exclude[i] << std::endl;
    }
    //print no_exclude
    std::cout << "No_exclude_Request : " << std::endl;
    for (size_t i = 0; i < no_exclude.size(); i++)
    {
        std::cout << no_exclude[i] << std::endl;
    }
    std::cout << "------------------" << std::endl;
    // print R arrival_city
    std::cout << "R Arrival_city : " << R["drop"][EVENT_CITY_INDEX] << std::endl;
    // print R1 arrival_city
    std::cout << "R1 Arrival_city : " << R1["drop"][EVENT_CITY_INDEX] << std::endl;

    return 0;
}

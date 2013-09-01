#include <iostream>
#include <fstream>
#include <stdexcept>
#include <regex>
#include <map>

#include <gtest/gtest.h>

#include <PunyHeaders/common_functions.hpp>

using namespace std;

size_t findDAColumn(const vector<string> & column_names) {
    const regex DA_column("DA:.*");
    size_t column_no = 0;
    for (const string & name : column_names) {
        if (regex_match(name, DA_column)) {
            return column_no;
        }
        column_no++;
    }

    throw runtime_error("DA_column not found.\n");
}

bool isComponent(const string & column_name) {
    const regex component_expr("((DV:.*)|(TR:.*(i|a))|(TR:.*:(Inhibitors|Activators)))");
    return regex_match(column_name, component_expr);
}

string obtainName(const string & column_name) {
    if (regex_match(column_name, regex("DV:.*")))
        return column_name.substr(3, column_name.size() - 3);
    if (regex_match(column_name, regex("TR:.*(i|a)")))
        return column_name.substr(3, column_name.size() - 4);
    if (regex_match(column_name, regex("TR:.*:Inhibitors")))
        return column_name.substr(3, column_name.size() - 3 - strlen("Inhibitors") - 1);
    if (regex_match(column_name, regex("TR:.*:Activators")))
        return column_name.substr(3, column_name.size() - 3 - strlen("Activators") - 1);
    return "";
}

TEST(FunctionTest, NameForming) {
    EXPECT_STREQ("Test", obtainName("DV:Test").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Testa").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Testi").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Test:Inhibitors").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Test:Activators").c_str());
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        ::testing::InitGoogleTest( &argc, argv );
        return RUN_ALL_TESTS();
    }

    fstream input_file(argv[1], ios::in);
    if (!input_file)
        throw invalid_argument("Wrong filename \"" + string(argv[1]) + "\".\n");

    string names_line;
    getline(input_file, names_line);

    vector<string> column_names;
    boost::split(column_names, names_line, boost::is_any_of(","));

    size_t DA_colum = findDAColumn(column_names);

    map<size_t, string> components;
    size_t column_no = 0;
    for (const string & column : column_names) {
        if (isComponent(column)) {
            components.insert(make_pair(column_no, obtainName(column)));
        }
        column_no++;
    }

    vector<vector<vector<double> > > timepoints; // Measurements are grouped by a timepoint, not by the experiment!
    timepoints.push_back(vector<vector<double> >());
    string measurement_time = "0";
    size_t time_point = 0;
    string data_line;
    while(getline(input_file, data_line)) {
        vector<string> data;
        boost::split(data, data_line, boost::is_any_of(","));
        if (measurement_time.compare(data[DA_colum]) != 0) {
            measurement_time = data[DA_colum];
            time_point++;
            timepoints.push_back(vector<vector<double> >());
        }
        vector<double> values;
        for (const pair<size_t, string> & component : components) {
            values.push_back(boost::lexical_cast<double>(data[component.first]));
        }
        timepoints[time_point].push_back(values);
    }

    return 0;
}


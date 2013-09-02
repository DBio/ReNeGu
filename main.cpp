#include <iostream>
#include <fstream>
#include <stdexcept>
#include <regex>
#include <map>
#include <cmath>

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

bool isComponent(const string & column_name) {
    return !obtainName(column_name).empty();
}

bool isMeasured(const string & column_name) {
    const regex component_expr("DV:.*");
    return regex_match(column_name, component_expr);
}

TEST(FunctionTest, NameForming) {
    EXPECT_TRUE(isComponent("DV:Test"));
    EXPECT_TRUE(isComponent("TR:Testa"));
    EXPECT_TRUE(isComponent("TR:Testi"));
    EXPECT_TRUE(isComponent("TR:Test:Inhibitors"));
    EXPECT_TRUE(isComponent("TR:Test:Activators"));
    EXPECT_FALSE(isComponent("VD:Test"));
    EXPECT_FALSE(isComponent("TR:TestA"));

    EXPECT_TRUE(isMeasured("DV:Test"));
    EXPECT_FALSE(isMeasured("VD:Test"));
    EXPECT_FALSE(isMeasured("TR:Testi"));

    EXPECT_STREQ("Test", obtainName("DV:Test").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Testa").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Testi").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Test:Inhibitors").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Test:Activators").c_str());
}

struct ComponentData {
    size_t column_no;
    string name;
    bool measured; ///< True iff the component is not either stimulated or inhibited.
    set<string> inhibitors; ///< Names of the species that inhibit this node.
    set<string> activators; ///< Names of the species that activate this node.
};

vector<ComponentData> getComponenets(const vector<string> & column_names) {
    vector<ComponentData>  components;
    size_t column_no = 0;
    for (const string & column : column_names) {
        if (isComponent(column)) {
            components.push_back({column_no, obtainName(column), isMeasured(column), set<string>(), set<string>()});
        }
        column_no++;
    }
    return components;
}

vector<vector<vector<double> > > getTimepoints(const vector<ComponentData> & components, const size_t DA_colum, fstream & input_file) {
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
        for (const ComponentData & component : components) {
            values.push_back(boost::lexical_cast<double>(data[component.column_no]));
        }
        timepoints[time_point].push_back(values);
    }
    return timepoints;
}

vector<vector<vector<double> > > normalize(const vector<vector<vector<double> > > & timepoints) {
    vector<vector<vector<double> > > normalized = timepoints;

    const size_t comp_count = normalized.front().front().size();
    vector<double> minimals(comp_count, numeric_limits<double>::max());
    vector<double> maximals(comp_count, 0.0);
    for (const vector<vector<double> > & timepoint : timepoints) {
        for (const vector<double> & measurement : timepoint) {
            for (const size_t comp_no : scope(measurement)) {
                minimals[comp_no] = min(minimals[comp_no], measurement[comp_no]);
                maximals[comp_no] = max(maximals[comp_no], measurement[comp_no]);
            }
        }
    }

    vector<double> difference;
    for (const size_t comp_no : range(comp_count)) {
        difference.push_back(maximals[comp_no] - minimals[comp_no]);
    }

    for (vector<vector<double> > & timepoint : normalized) {
        for (vector<double> & measurement : timepoint) {
            for (const size_t comp_no : scope(measurement)) {
                measurement[comp_no] = (measurement[comp_no] - minimals[comp_no]) / difference[comp_no];
            }
        }
    }

    return normalized;
}

class DataTest : public ::testing::Test {
protected:
    vector<vector<vector<double> > > data1;

    void SetUp() override {
        data1 = {{{0,0.4},{1,0.6}},{{0.5,0.8}, {0.25,0.5}}};
    }
};

TEST_F(DataTest, Normalize) {
    vector<vector<vector<double> > > values = normalize(data1);
    EXPECT_DOUBLE_EQ(0., values[0][0][0]);
    EXPECT_DOUBLE_EQ(0., values[0][0][1]);
    EXPECT_DOUBLE_EQ(1., values[0][1][0]);
    EXPECT_DOUBLE_EQ(0.5, values[0][1][1]);
    EXPECT_DOUBLE_EQ(0.5, values[1][0][0]);
    EXPECT_DOUBLE_EQ(1., values[1][0][1]);
    EXPECT_DOUBLE_EQ(0.25, values[1][1][0]);
    EXPECT_DOUBLE_EQ(0.25, values[1][1][1]);
}

struct Interval {
    const vector<double> & first_;
    const vector<double> & second_;
    vector<double> changes;
    vector<double> productions;

    NO_COPY_SHORT(Interval)
    DEFAULT_MOVE(Interval)

    Interval(const vector<double> & first, const vector<double> & second) : first_(first), second_(second) {
        double avg_production = 0;
        size_t production_cases = 0;

        // Compute change
        for (const size_t val_no : scope(first)) {
            const double change = second[val_no] - first[val_no];
            changes.push_back(change);
            if (step > 0.) {
                avg_production += change;
                production_cases++;
            }
        }
        avg_production /= avg_production;

        // Compute production factor - positive if bigger than average, negative otherwise
        for (const double change : changes) {
            if (change <= 0.)
                productions.push_back(0.);
            else
                productions.push_back(change - avg_production);
        }
    }
};

vector<Interval> makeIntervals(const vector<vector<vector<double> > > & values) {
    vector<Interval> intervals;

    for (const size_t experiment_no : scope(values[0])) {
        for (const size_t measurement_no : range(values.size() - 1)) {
            intervals.push_back(Interval(values[measurement_no][experiment_no], values[measurement_no+1][experiment_no]));
        }
    }

    return intervals;
}

TEST_F(DataTest, MakeIntervals) {
    vector<Interval> intervals = makeIntervals(normalize(data1));
    EXPECT_DOUBLE_EQ(0.5, intervals[0].change[0]);
    EXPECT_DOUBLE_EQ(1., intervals[0].change[1]);
    EXPECT_DOUBLE_EQ(-0.75, intervals[1].change[0]);
    EXPECT_DOUBLE_EQ(-0.25, intervals[1].change[1]);
}

double inline logistic(const double x) {
    return 1 / (1 + exp(-x));
}

double inline logistic(const double x, const double y) {
    return logistic(x) + logistic(y) - logistic(x)*logistic(y);
}

double inline logistic(const double x, const double y, const double z) {
    return logistic(x) + logistic(y) + logistic(z)
           - logistic(x)*logistic(y) - logistic(x)*logistic(z) - logistic(y)*logistic(z)
           + logistic(x)*logistic(y)*logistic(z);
}


vector<ComponentData> computeRegulators(const vector<ComponentData> & components, vector<Interval> & intervals) {
    vector<ComponentData> regulated = components;

    for (const size_t target_no : scope(components)) {
        map<double, set<string> > mutual_inf;
        if (!components[target_no].measured)
            continue;

        // Test single inhibitors.
        for (const size_t source_no : scope(components)) {
            double effect = 0.;
            for (const Interval & interval : intervals) {
                effect += -1 * interval.changes[target_no] * logistic(interval.first_[source_no]);
            }
            effect /= intervals.size();
            mutual_inf.insert(make_pair(effect, set<string>({components[source_no].name})));
        }

        // Test two inhibitors.
        for (const size_t source_1_no : scope(components)) {
            for (const size_t source_2_no : range(source_1_no + 1, components.size())) {
                double effect = 0.;
                for (const Interval & interval : intervals) {
                    effect += -1 * interval.changes[target_no] * logistic(interval.first_[source_1_no],interval.first_[source_2_no]);
                }
                effect /= intervals.size();
                mutual_inf.insert(make_pair(effect, set<string>({components[source_1_no].name, components[source_2_no].name})));
            }
        }

        // Test three inhibitors.
        for (const size_t source_1_no : scope(components)) {
            for (const size_t source_2_no : range(source_1_no + 1, components.size())) {
                for (const size_t source_3_no : range(source_1_no + 1, components.size())) {
                    double effect = 0.;
                    for (const Interval & interval : intervals) {
                        effect += -1 * interval.changes[target_no] * logistic(interval.first_[source_1_no],interval.first_[source_2_no],interval.first_[source_3_no]);
                    }
                    effect /= intervals.size();
                    mutual_inf.insert(make_pair(effect, set<string>({components[source_1_no].name, components[source_2_no].name, components[source_3_no].name})));
                }
            }
        }


        auto max_it = max_element(mutual_inf.begin(), mutual_inf.end());
        // Pick the one with the biggest positive influence - if there is none, just take empty.
        regulated[target_no].inhibitors = max_it->first > 0 ? max_it->second : set<string>();
    }

    return regulated;
}

int main(int argc, char* argv[]) {
    // if (argc < 2) {
        ::testing::InitGoogleTest( &argc, argv );
        RUN_ALL_TESTS();
    // }

    fstream input_file(argv[1], ios::in);
    if (!input_file)
        throw invalid_argument("Wrong filename \"" + string(argv[1]) + "\".\n");

    string names_line;
    getline(input_file, names_line);

    vector<string> column_names;
    boost::split(column_names, names_line, boost::is_any_of(","));

    size_t DA_colum = findDAColumn(column_names);
    vector<ComponentData> components = getComponenets(column_names);
    vector<vector<vector<double> > > timepoints = getTimepoints(components, DA_colum, input_file);
    timepoints = normalize(timepoints);
    vector<Interval> intervals = makeIntervals(timepoints);
    components = computeRegulators(components, intervals);

    return 0;
}


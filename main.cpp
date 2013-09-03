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
    if (regex_match(column_name, regex("TR:.*:Inhibitors")))
        return column_name.substr(3, column_name.size() - 3 - strlen("Inhibitors") - 1);
    if (regex_match(column_name, regex("TR:.*:Stimuli")))
        return column_name.substr(3, column_name.size() - 3 - strlen("Stimuli") - 1);
    if (regex_match(column_name, regex("TR:.*(i|a)")))
        return column_name.substr(3, column_name.size() - 4);
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
    EXPECT_TRUE(isComponent("TR:Test:Stimuli"));
    EXPECT_FALSE(isComponent("VD:Test"));
    EXPECT_FALSE(isComponent("TR:TestA"));

    EXPECT_TRUE(isMeasured("DV:Test"));
    EXPECT_FALSE(isMeasured("VD:Test"));
    EXPECT_FALSE(isMeasured("TR:Testi"));

    EXPECT_STREQ("Test", obtainName("DV:Test").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Testa").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Testi").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Test:Inhibitors").c_str());
    EXPECT_STREQ("Test", obtainName("TR:Test:Stimuli").c_str());
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

vector<vector<string> > getData(fstream & input_file) {
    vector<vector<string > > data;

    string line;
    while(getline(input_file, line)) {
        vector<string> line_data;
        boost::split(line_data, line, boost::is_any_of(","));
        data.push_back(line_data);
    }

    return data;
}

inline size_t getTime(const vector<string> & line, const size_t DA_column) {
    return boost::lexical_cast<size_t>(line[DA_column]);
}

set<size_t> getTimes(const vector<vector<string> > & data, const size_t DA_column) {
    set<size_t> times;

    for (const vector<string> line : data) {
        times.insert(getTime(line, DA_column));
    }

    return times;
}

inline bool isInTime(const vector<string> line, const size_t DA_column, const size_t time) {
    return (getTime(line, DA_column) == time);
}

bool matchesSetting(const vector<ComponentData> & components, const vector<string> & source, const vector<string> target) {
    for (const size_t comp_no : scope(components)) {
        if (components[comp_no].measured)
            continue;

        if (source[comp_no].compare(source[comp_no]) != 0)
            return false;
    }

    return true;
}

vector<pair<vector<string>, vector<string> > > getTimesteps(const vector<ComponentData> & components, const size_t DA_column, const vector<vector<string> > & data) {
    vector<pair<vector<string>, vector<string> > > timesteps; // Measurements are grouped by a timepoint, not by the experiment!

    const set<size_t> times = getTimes(data, DA_column);

    for (const vector<string> & source_line : data) {
        const size_t source_time = getTime(source_line, DA_column);
        auto time_it = times.find(source_time);
        if(++time_it == times.end())
            continue;

        for (const vector<string> & target_line : data) {
            if (!isInTime(target_line, DA_column, *time_it))
                continue;

            if (matchesSetting(components, source_line, target_line)) {
                timesteps.push_back(make_pair(source_line, target_line));
                break;
            }
        }
    }

    return timesteps;
}


vector<pair<vector<double>, vector<double> > > convert(const vector<ComponentData> & components, const vector<pair<vector<string>, vector<string> > > & timesteps) {
    vector<pair<vector<double>, vector<double> > > converted;

    for (const pair<vector<string>, vector<string> >  & timestep : timesteps) {
        vector<double> first;
        vector<double> second;
        for (const ComponentData & component : components) {
            first.push_back(boost::lexical_cast<double>(timestep.first[component.column_no]));
            second.push_back(boost::lexical_cast<double>(timestep.second[component.column_no]));
        }
        converted.push_back(make_pair(first, second));
    }

    return converted;
}

vector<pair<vector<double>, vector<double> > > normalize(const vector<pair<vector<double>, vector<double> > > & timesteps) {
    vector<pair<vector<double>, vector<double> > > normalized;

    if (timesteps.empty())
        throw runtime_error("Nothing to normalize.");
    const size_t comp_count = timesteps.front().first.size();
    vector<double> minimals(comp_count, numeric_limits<double>::max());
    vector<double> maximals(comp_count, 0.0);
    for (const pair<vector<double>, vector<double> >  & timestep : timesteps) {
        for (const size_t comp_no : range(comp_count)) {
            minimals[comp_no] = min(minimals[comp_no], min(timestep.first[comp_no], timestep.second[comp_no]));
            maximals[comp_no] = max(maximals[comp_no], max(timestep.first[comp_no], timestep.second[comp_no]));
        }
    }

    vector<double> difference;
    for (const size_t comp_no : range(comp_count)) {
        difference.push_back(maximals[comp_no] - minimals[comp_no]);
    }

    for (const pair<vector<double>, vector<double> >  & timestep : timesteps) {
        vector<double> first, second;
        for (const size_t comp_no : range(comp_count)) {
            first.push_back((timestep.first[comp_no] - minimals[comp_no]) / difference[comp_no]);
            second.push_back((timestep.second[comp_no] - minimals[comp_no]) / difference[comp_no]);
        }
        normalized.push_back(make_pair(first, second));
    }

    return normalized;
}

class DataTest : public ::testing::Test {
protected:
    vector<pair<vector<double>, vector<double> > > data1;

    void SetUp() override {
        data1 = {{{0,0.4},{0.5,0.8}},{{1,0.6}, {0.25,0.5}}};
    }
};

//TEST_F(DataTest, Normalize) {
//    vector<vector<vector<double> > > values = normalize(data1);
//    EXPECT_DOUBLE_EQ(0., values[0][0][0]);
//    EXPECT_DOUBLE_EQ(0., values[0][0][1]);
//    EXPECT_DOUBLE_EQ(1., values[0][1][0]);
//    EXPECT_DOUBLE_EQ(0.5, values[0][1][1]);
//    EXPECT_DOUBLE_EQ(0.5, values[1][0][0]);
//    EXPECT_DOUBLE_EQ(1., values[1][0][1]);
//    EXPECT_DOUBLE_EQ(0.25, values[1][1][0]);
//    EXPECT_DOUBLE_EQ(0.25, values[1][1][1]);
//}

struct Interval {
    const vector<double> first_;
    const vector<double> second_;
    vector<double> changes_;
    vector<double> productions_;

    NO_COPY_SHORT(Interval)
    DEFAULT_MOVE(Interval)

    Interval(const vector<double> & first, const vector<double> & second) : first_(first), second_(second) {
        double avg_production = 0;
        size_t production_cases = 0;

        // Compute change
        for (const size_t val_no : scope(first)) {
            const double change = second[val_no] - first[val_no];
            changes_.push_back(change);
            if ((change) > 0.) {
                avg_production += change;
                production_cases++;
            }
        }
        avg_production /= production_cases;

        // Compute production factor - positive if bigger than average, negative otherwise
        for (const double change : changes_) {
            if (change <= 0.)
                productions_.push_back(0.);
            else
                productions_.push_back(change - avg_production);
        }
    }
};

vector<Interval> makeIntervals(const vector<pair<vector<double>, vector<double> > > & timesteps) {
    vector<Interval> intervals;

    for (const pair<vector<double>, vector<double> >  & timestep : timesteps) {
        intervals.push_back(Interval(timestep.first, timestep.second));
    }

    return intervals;
}

TEST_F(DataTest, MakeIntervals) {
    vector<Interval> intervals = makeIntervals(normalize(data1));
    EXPECT_DOUBLE_EQ(0.5, intervals[0].changes_[0]);
    EXPECT_DOUBLE_EQ(1., intervals[0].changes_[1]);
    EXPECT_DOUBLE_EQ(-0.75, intervals[1].changes_[0]);
    EXPECT_DOUBLE_EQ(-0.25, intervals[1].changes_[1]);
}

double inline logistic(const double x) {
    return 1 / (1 + exp(-(x*12-6)));
}

double inline logistic(const double x, const double y) {
    return logistic((x+y)/2);
}

double inline logistic(const double x, const double y, const double z) {
    return logistic((x+y+z)/3);
}

TEST(FunctionTest, logistic) {
    cout << "logistic(0.5) " << logistic(0.5) << ", logistic(0.5, 0.0) " << logistic(0.5, 0.0) << ", logistic(0.5, 0.0, 0.0) " << logistic(0.5, 0.0, 0.0) << ".\n";
    cout << "logistic(0.5) " << logistic(0.5) << ", logistic(0.5, 0.5) " << logistic(0.5, 0.5) << ", logistic(0.5, 0.5, 0.5) " << logistic(0.5, 0.5, 0.5) << ".\n";
    cout << "logistic(0.5) " << logistic(0.5) << ", logistic(0.5, 1.0) " << logistic(0.5, 1.0) << ", logistic(0.5, 1.0, 1.0) " << logistic(0.5, 1.0, 1.0) << ".\n";
    cout << "logistic(1.0) " << logistic(1.0) << ", logistic(1.0, 1.0) " << logistic(1.0, 1.0) << ", logistic(1.0, 1.0, 1.0) " << logistic(1.0, 1.0, 1.0) << ".\n";
    cout << "logistic(0.1,0.9) " << logistic(0.1,0.9) << "logistic(0.5,0.5) "  << logistic(0.5,0.5) << ".\n";
}

vector<ComponentData> computeRegulators(const vector<ComponentData> & components, vector<Interval> & intervals, bool inhibition) {
    vector<ComponentData> regulated = components;
    const size_t total = pow(components.size(), 3) * count_if(components.begin(), components.end(), [](const ComponentData & comp){return comp.measured;}) * 2;
    size_t current = inhibition ? 0 : total / 2;

    for (const size_t target_no : scope(components)) {
        map<set<string>,  double> mutual_inf;
        if (!components[target_no].measured)
            continue;

        // Test all combinations up to three regulators (combination of multiple selfs == just single one).
        for (const size_t source_1_no : scope(components)) {
            for (const size_t source_2_no : scope(components)) {
                for (const size_t source_3_no : scope(components)) {
                     cout << "Testing: " << current++ << "/" << total << "     \r";
                    double effect = 0.;
                    for (const Interval & interval : intervals) {
                        double factor = inhibition ? -1 * interval.changes_[target_no] : interval.productions_[target_no];
                        effect += factor * logistic(interval.first_[source_1_no],interval.first_[source_2_no],interval.first_[source_3_no]);
                    }
                    mutual_inf.insert(make_pair(set<string>({components[source_1_no].name, components[source_2_no].name, components[source_3_no].name}), effect));
                }
            }
        }

        auto max_it = max_element(mutual_inf.begin(), mutual_inf.end(),
                                  [](const pair<set<string>, double> & A, const pair<set<string>, double> & B){
            return A.second < B.second;
        });
        // Pick the one with the biggest positive influence - if there is none, just take empty.

        set<string> & sources = inhibition ? regulated[target_no].inhibitors : regulated[target_no].activators;
        sources = (max_it->second - 0.05 * max_it->first.size()) > 0 ? max_it->first : set<string>();
    }

    return regulated;
}

string getModelName(const string & argument) {
    string name;
    for (const char ch : argument) {
        if (ch == '.')
            break;
        name.push_back(ch);
        if (ch == '/' || ch == '\\')
            name.clear();
    }
    return name;
}


void output(const vector<ComponentData> & components, const string & file_name) {
    ofstream fout(file_name + ".sif", ios::out);

    for (const ComponentData & component: components) {
        for (const string & source : component.activators) {
            fout << source << " 1 " << component.name << endl;
        }
        for (const string & source : component.inhibitors) {
            fout << source << " -1 " << component.name << endl;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        ::testing::InitGoogleTest( &argc, argv );
        return RUN_ALL_TESTS();
    }

    string filename(argv[1]);
    fstream input_file(filename, ios::in);
    if (!input_file)
        throw invalid_argument("Wrong filename \"" + filename + "\".\n");

    string names_line;
    getline(input_file, names_line);

    vector<string> column_names;
    boost::split(column_names, names_line, boost::is_any_of(","));

    size_t DA_colum = findDAColumn(column_names);
    vector<ComponentData> components = getComponenets(column_names);
    vector<vector<string> > data = getData(input_file);
    vector<pair<vector<string>, vector<string> > > timepoints = getTimesteps(components, DA_colum, data);
    vector<pair<vector<double>, vector<double> > > converted = convert(components, timepoints);
    vector<pair<vector<double>, vector<double> > > normalized = normalize(converted);
    vector<Interval> intervals = makeIntervals(normalized);
    components = computeRegulators(components, intervals, true);
    components = computeRegulators(components, intervals, false);

    output(components, getModelName(filename));

    return 0;
}


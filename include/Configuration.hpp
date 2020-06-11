#ifndef __CONFIGURATION__
#define __CONFIGURATION__

#include <map>
#include <string>
//#include <fstream>
//#include <iostream>

using namespace std;

/// This class can read and keep values of any configuration file written in a format like this:
/// #
/// # Example parameters
/// #
///
/// int_param     = 2 # NOTE: Can write comment after param
/// float_param   = 25.001
/// sci_param     = 1.013e-14
/// bool_param    = false
/// string_param  = a string
///
/// Example usage:
///
///Configuration   config;
///int             int_param;
///double          float_param;
///double          sci_param;
///bool            bool_param;
///string          string_param;
/// 
///config.Load("./simulation.cfg");
/// 
///if (config.Get("int_param", int_param)        &&
///    config.Get("float_param", float_param)    &&
///    config.Get("sci_param",    sci_param)     &&
///    config.Get("bool_param",   bool_param)    &&
///    config.Get("string_param", string_param))
///{
///      ...
///}
///else
///{
///  cout << "Missing parameter in configuration file." << endl;
///}

class Configuration
{
  public:
    // clear all values
    void Clear();

    // load a configuration file
    bool Load(const string& File);

    // check if value associated with given key exists
    bool Contains(const string& key) const;

    // get value associated with given key
    bool Get(const string& key, string& value) const;
    bool Get(const string& key, int&    value) const;
    bool Get(const string& key, long&   value) const;
    bool Get(const string& key, double& value) const;
    bool Get(const string& key, bool&   value) const;

  private:
    // the container
    map<string,string> data;

    // remove leading and trailing tabs and spaces
    static string Trim(const string& str);
};

#endif // __CONFIGURATION__

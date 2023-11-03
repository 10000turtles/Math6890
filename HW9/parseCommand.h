#ifndef PARSE_COMMAND_H
#define PARSE_COMMAND_H

#include <string>
using std::string;

//-----------------------------------------------------------------------------------
// Parseacommandfromastring"line"andreturn ainteger"value"
// Example:
// line="-Nx=10"(input)
// command="-Nx="(input)
// value=(output)
// Optionalinput:echo:echo=truemeansprintwhenavariableischanged.
// Returnvalue:1ifcommandwasfoundinline
// 0ifnotfound
//-----------------------------------------------------------------------------------
int parseCommand(const string &line, const string &command, int &value, bool echo = true)
{
    int len = command.length();
    if (line.substr(0, len) == command)
    {
        sscanf(line.substr(len).c_str(), "%d", &value);
        if (echo)
            printf("parseCommand:SETTING%s%d\n", command.c_str(), value);
        return 1;
    }
    return 0;
}

// Parseadouble
int parseCommand(const string &line, const string &command, double &value, bool echo = true)
{
    int len = command.length();
    if (line.substr(0, len) == command)
    {
        sscanf(line.substr(len).c_str(), "%le", &value);
        if (echo)
            printf("parseCommand:SETTING%s%e\n", command.c_str(), value);
        return 1;
    }
    return 0;
}

// Parseastring
int parseCommand(const string &line, const string &command, string &value, bool echo = true)
{
    int len = command.length();
    if (line.substr(0, len) == command)
    {
        value = line.substr(len);
        if (echo)
            printf("parseCommand:SETTING%s%s\n", command.c_str(), value.c_str());
        return 1;
    }
    return 0;
}

#endif
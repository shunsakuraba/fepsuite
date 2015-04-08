#include <string>

inline
std::string rtrim(const std::string& s)
{
  size_t p = s.find_last_not_of(' ');
  size_t len = p + 1;
  return s.substr(0, len);
}

inline
std::string ltrim(const std::string& s)
{
  size_t p = s.find_first_not_of(' ');
  return s.substr(p, std::string::npos);
}

inline
std::string trim(const std::string& s)
{
  return ltrim(rtrim(s));
}

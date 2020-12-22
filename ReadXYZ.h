#ifndef READXYZ_H
#define READXYZ_H

#include <string>
#include <fstream>
#include <sstream>
#include <cctype>
#include <vector>

template <typename Point3>
inline bool ReadXYZ(const std::string& name,
        std::vector<Point3>& points, std::vector<Point3>& normals)
{
    std::ifstream in(name);
    if(!in) return false;

    points.clear();

    std::string line;
    std::stringstream ss;
    Point3 p;
    Point3 n;

    bool first = true;
    bool has_normal = false;
    while(!in.eof())
    {
        std::getline(in, line);
        if(in.bad()) return false;

        ss.str(line);
        ss.clear();
        ss >> std::ws;
        if(ss.eof()) continue;
        auto c = ss.peek();
        if(!std::isdigit(c) && c != '-' && c != '+') continue;

        auto getPoint = [&](Point3& p) -> bool {
            for(int i = 0; i < 3; i++)
            {
                if(ss.eof()) return false;
                ss >> p[i];
            }
            return true;
        };

        if(!getPoint(p)) continue;

        if(first)
        {
            first = false;
            if(!ss.eof())
            {
				ss >> std::ws;
                c = ss.peek();
                if(std::isdigit(c) || c == '-' || c == '+')
                    has_normal = true;
            }
        }

        if(has_normal)
        {
            if(!getPoint(n)) continue;
        }

        points.emplace_back(p);
        if(has_normal)
            normals.emplace_back(n);
    }
    return true;
}
#endif

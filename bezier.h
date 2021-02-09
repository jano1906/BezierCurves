#ifndef JO418361_BEZIER_H
#define JO418361_BEZIER_H

#include <iostream>
#include <utility>
#include <vector>
#include <functional>
#include <cassert>
#include <cmath>

namespace bezier::types {
    using real_t = double;
    using node_index_t = uint32_t;
    class point_2d {
    public:
        const real_t X;
        const real_t Y;
        point_2d(real_t x, real_t y) : X(x), Y(y){}
    };

    point_2d operator+(point_2d p1, point_2d p2){
        return {p1.X + p2.X, p1.Y + p2.Y};
    }
    point_2d operator*(point_2d p, real_t c){
        return {p.X*c, p.Y*c};
    }
    point_2d operator*(real_t c, point_2d p){
        return p*c;
    }
    bool operator==(point_2d p1, point_2d p2){
        return p1.X==p2.X && p1.Y==p2.Y;
    }
}

std::ostream& operator<<(std::ostream& os, bezier::types::point_2d p){
    return os<<"("<<p.X<<", "<<p.Y<<")";
}

namespace bezier::constants {
    uint16_t constexpr NUM_OF_CUBIC_BEZIER_NODES = 4;
    const double ARC = 4*(std::sqrt(2.) - 1)/3;
}

namespace {
    template<typename ID, typename T>
    std::function<T (ID)> fun(ID n, T x){
        return [x, n](ID id) {
                if(id == n)
                    return x;
                throw std::out_of_range("a curve node index is out of range");
            };
    }

    template<typename ID, typename T, typename ...Args>
    std::function<T (ID)> fun(ID n, T x, Args... args) {
        return [n, x, args...](ID id) {
            if(id == n)
                return x;
            if(id > n)
                return fun(n+1, args...)(id);
            throw std::out_of_range("a curve node index is out of range");
        };
    }

    using namespace bezier::types;
    using Matrix = real_t[2][2];
    point_2d matmul(Matrix m, point_2d p){
        return {p.X*m[0][0] + p.Y*m[0][1], p.X*m[1][0] + p.Y*m[1][1]};
    }
}
namespace bezier {
    using namespace bezier::types;
    using namespace bezier::constants;
    using curve_t = std::function<point_2d (node_index_t)>;

    curve_t Cup() {
        point_2d p[4]{{-1,1}, {-1,-1}, {1,-1}, {1,1}};
        return fun(0, p[0],p[1],p[2],p[3]);
    }
    curve_t Cap() {
        point_2d p[4]{{-1,-1}, {-1,1}, {1,1}, {1,-1}};
        return fun(0, p[0],p[1],p[2],p[3]);
    }
    curve_t ConvexArc() {
        point_2d p[4]{{0,1}, {ARC,1}, {1,ARC}, {1,0}};
        return fun(0, p[0],p[1],p[2],p[3]);
    }
    curve_t ConcaveArc() {
        point_2d p[4]{{0,1}, {0,1-ARC}, {1-ARC,0}, {1,0}};
        return fun(0, p[0],p[1],p[2],p[3]);
    }
    curve_t LineSegment(point_2d p, point_2d q) {
        return fun(0, p,p,q,q);
    }
    curve_t MovePoint(const curve_t& f, node_index_t i, real_t x, real_t y) {
        return [f,i,x,y](node_index_t id){
              if(id == i) {
                  point_2d p{x,y};
                  return f(id) + p;
              }
            return f(id);
        };
    }
    curve_t Rotate(const curve_t& f, real_t a) {
        //converting to radians
        a = a * M_PI / 180;
        return [f, a] (node_index_t id){
            Matrix rot = {{cos(a), -sin(a)},
                          {sin(a), cos(a)}};
            return matmul(rot, f(id));
        };
    }
    curve_t Scale(const curve_t& f, real_t xsc, real_t ysc) {
        return [f, xsc, ysc] (node_index_t id){
            Matrix sca = {{xsc, 0},
                          {0, ysc}};
            return matmul(sca, f(id));
        };
    }
    curve_t Translate(const curve_t& f, real_t xtr, real_t ytr) {
        point_2d tr = {xtr, ytr};
        return [f, tr] (node_index_t id){
            return f(id) + tr;
        };
    }
    curve_t Concatenate(const curve_t& c1, const curve_t& c2) {
        return[c1, c2](node_index_t id) {
            if(id >= 2*NUM_OF_CUBIC_BEZIER_NODES)
                throw std::out_of_range("a curve node index is out of range");
            curve_t c[2] = {c1,c2};
            return c[id/NUM_OF_CUBIC_BEZIER_NODES](id%NUM_OF_CUBIC_BEZIER_NODES);
        };
    }
    template<typename... Args>
    curve_t Concatenate(const curve_t& c, Args... args) {
        return[c, args...] (node_index_t id) {
            if(id < NUM_OF_CUBIC_BEZIER_NODES)
                return c(id);
            return Concatenate(args...)(id - NUM_OF_CUBIC_BEZIER_NODES);
        };
    }
    class P3CurvePlotter {
        using res_t = uint32_t;
        using seg_t = uint32_t;
        using plot_t = const std::vector<std::vector<bool>>;
        res_t res;
        seg_t seg;
        curve_t curve;
        real_t sq_x_len;
        real_t sq_y_len;
        plot_t plot;

        bool inSquare(point_2d p) const {
            return  -sq_x_len/2<=p.X && p.X<=sq_x_len/2 &&
                    -sq_y_len/2<=p.Y && p.Y<=sq_y_len/2;
        }

        std::pair<res_t, res_t> getPixelCoords(point_2d p) const {
            return {(res-1)*(p.X/sq_x_len + 1./2), (res-1)*(p.Y/sq_y_len + 1./2)};
        }

        std::vector<std::vector<bool>> computePlot() {
            std::vector<std::vector<bool>> plt(res, std::vector<bool>(res, false));
            for(res_t x=0; x <= res*res; x++) {
                for(seg_t s = 0; s<seg; s++) {
                    real_t t = ((real_t)x)/(res*res);
                    point_2d p = operator()(curve, t, s);
                    if(inSquare(p)) {
                        auto coords = getPixelCoords(p);
                        plt[coords.first][coords.second] = true;
                    }
                }
            }
            return plt;
        }
    public:
        explicit P3CurvePlotter(curve_t curve, seg_t seg = 1, res_t res = 80,
                                real_t sq_x_len = 2., real_t sq_y_len = 2.)
            : curve(std::move(curve)), seg(seg), res(res),
                sq_x_len(sq_x_len), sq_y_len(sq_y_len), plot(computePlot()) {}

        void Print(std::ostream& os = std::cout, char fb = '*', char bg = ' ') const {
            for(res_t y = 0; y<res; y++) {
                for(res_t x = 0; x<res;x++) {
                    if(plot[x][res-1-y]) {
                        os<<fb;
                    } else {
                        os<<bg;
                    }
                }
                os<<std::endl;
            }
        }
        point_2d operator()(const curve_t& c, real_t t, seg_t s) const {
            assert(0<=t && t<=1);

            node_index_t start = s*NUM_OF_CUBIC_BEZIER_NODES;
            std::vector<point_2d> b;
            for(node_index_t i = start; i < start+NUM_OF_CUBIC_BEZIER_NODES;i ++)
                b.push_back(c(i));

            for(int k = NUM_OF_CUBIC_BEZIER_NODES; k>0; k--){
                for(int i=1;i<k;i++){
                    b.push_back(b[b.size()-k]*(1-t) + b[b.size()-k+1]*t);
                }
            }
            return b[b.size()-1];
        }
    };
}

#endif //JO418361_BEZIER_H

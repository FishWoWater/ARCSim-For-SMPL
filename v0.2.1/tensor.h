#ifndef TENSOR_HPP
#define TENSOR_HPP

#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace std;

template<int T>
class TensorD : public Eigen::Tensor<float, T>{
private:
    int tDim = -1;
public:
    typedef Eigen::Tensor<float, T> Base;
    typedef Eigen::TensorMap<Eigen::Tensor<float, T>> Map;
    typedef Eigen::TensorMap<Eigen::Tensor<float, 2>> Map2D;
    typedef Eigen::TensorMap<Eigen::Tensor<float, 3>> Map3D;
    typedef Eigen::TensorCwiseBinaryOp<Eigen::internal::scalar_sum_op<float>, const Base, const Base> BinaryOp;
    typedef Eigen::TensorContractionOp<const std::array<Eigen::IndexPair<int>,1ul>, const Base, const Eigen::Tensor<float, 2>> Contraction2Op;
    typedef Eigen::TensorContractionOp<const std::array<Eigen::IndexPair<int>,1ul>, const Base, const Map2D> Contraction2MapOp;

    TensorD(){
        tDim = T;
    }

    TensorD(const std::array<int, T>& size);

    TensorD( const Base &d );

    TensorD(const BinaryOp& x);

    TensorD(const Contraction2Op& x);

    TensorD(const Contraction2MapOp& x);

    TensorD<T> operator+(TensorD<T>& n);

    ~TensorD();

    void resize(const std::array<int, T>& size);

    void printSize();

    int depth();

    int rows();

    int cols();

    void printAtDepth(int d, int rowCount, int colCount);

    Eigen::MatrixXf getMatrix(int d);

    Eigen::MatrixXf getMatrix();

    void setFromMatrix(const Eigen::MatrixXf& m);

    TensorD<3> dot(const TensorD<2>& x);

    TensorD<3> dot(Eigen::MatrixXf& m);

    void print2D(int rowCount = 4, int colCount = 4);
    void print3D(int depthCount = 4, int rowCount = 4, int colCount = 4);
    void print();
};


std::string getSpaces(int spaces, bool empt);
void printEigen(const Eigen::MatrixXf& m, int colCount, int rowCount, int spaces);

#endif
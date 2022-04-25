#include "tensor.h"

template<int T>
TensorD<T>::TensorD(const std::array<int, T>& size){
    tDim = size.size();
    resize(size);
}

template<int T>
TensorD<T>::TensorD( const Base &d )
{
    tDim = d.dimensions().size();
}

template<int T>
TensorD<T>::TensorD(const BinaryOp& x)
{
    this->base() = x;
    tDim = this->base().dimensions().size();
}

template<int T>
TensorD<T>::TensorD(const TensorD<T>::Contraction2Op& x)
{
    this->base() = x;
    tDim = this->base().dimensions().size();
}

template<int T>
TensorD<T>::TensorD(const TensorD<T>::Contraction2MapOp& x)
{
    this->base() = x;
    tDim = this->base().dimensions().size();
}

template<int T>
TensorD<T> TensorD<T>::operator+(TensorD<T>& n){
    return (this->base()+n.base());
}

template<> TensorD<3>::~TensorD(){

}

template<int T>
TensorD<T>::~TensorD(){

}

template<int T>
void TensorD<T>::resize(const std::array<int, T>& size){
    tDim = size.size();
    this->Base::resize(size);
}

//template<int T>
template<> void TensorD<3>::resize(const std::array<int, 3>& size){
    tDim = size.size();
    this->Base::resize(size);
}

template<int T>
void TensorD<T>::printSize(){
        if(tDim == 3)
            cout << "D: " << this->dimensions()[0] << " R: " << this->dimensions()[1] << " C: " << this->dimensions()[2] << endl;
        else if(tDim == 2)
            cout << "R: " << this->dimensions()[0] << " C: " << this->dimensions()[1] << endl;
        else
            throw std::runtime_error("printSize dim > 3 not implemented. Or failed to call resize");
    }

template<int T>
int TensorD<T>::depth(){
    if(tDim != 3) throw std::runtime_error("This matrix has no depth");
    return this->dimensions()[0];
}

template<int T>
int TensorD<T>::rows(){
    if(tDim == 3)
        return this->dimensions()[1];
    else if(tDim == 2)
        return this->dimensions()[0];
}

template<int T>
int TensorD<T>::cols(){
    if(tDim == 3)
        return this->dimensions()[2];
    else if(tDim == 2)
        return this->dimensions()[1];
}

template<int T>
void TensorD<T>::printAtDepth(int d, int rowCount, int colCount){
    if(tDim != 3) throw std::runtime_error("This matrix has no depth");
    std::cout << std::setprecision(3) << std::fixed;
    cout << "[";
    for(int i=0; i<rows(); i++){
        if(i < rowCount || i >= rows()-rowCount){
            cout << "[";
            for(int j=0; j<cols(); j++){
                if(j < colCount || j >= cols()-colCount)
                    cout << this->operator()(d,i,j) << "     ";
                else if(j == colCount)
                    cout << " ... ";
            }
            cout << "]" << endl;
        }else if(i == rowCount){
            cout << "..." << endl;
            cout << "..." << endl;
        }
    }
    cout << "]" << endl;
}

template<int T>
Eigen::MatrixXf TensorD<T>::getMatrix(int d){
    Eigen::MatrixXf m;
    m.resize(rows(),cols());
    for(int r=0;r<rows(); r++){
        for(int c=0;c<cols(); c++){
            m(r,c) = this->operator()(d,r,c);
        }
    }
    return m;
}

template<int T>
Eigen::MatrixXf TensorD<T>::getMatrix(){
    Eigen::MatrixXf m;
    m.resize(rows(),cols());
    for(int r=0;r<rows(); r++){
        for(int c=0;c<cols(); c++){
            m(r,c) = this->operator()(r,c);
        }
    }
    return m;
}

template<int T>
void TensorD<T>::setFromMatrix(const Eigen::MatrixXf& m){
    if(!(m.rows() == rows() && m.cols() == cols())) resize({m.rows(), m.cols()});
    for(int r=0;r<rows(); r++){
        for(int c=0;c<cols(); c++){
            this->operator()(r,c) = m(r,c);
        }
    }
}

//template<int T>
template<> TensorD<3> TensorD<3>::dot(const TensorD<2>& x){
    Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>(2, 0) };
    return this->base().contract(x.base(), product_dims);
}

template<> TensorD<3> TensorD<3>::dot(Eigen::MatrixXf& m){
    auto map = Map2D(m.data(), m.rows(), m.cols());
    Eigen::array<Eigen::IndexPair<int>, 1> product_dims = { Eigen::IndexPair<int>(2, 0) };
    return this->base().contract(map, product_dims);
}


template<> void TensorD<2>::print2D(int rowCount, int colCount){
    printSize();
    std::cout << std::setprecision(3) << std::fixed;
    cout << "[";
    for(int i=0; i<rows(); i++){
        if(i < rowCount || i >= rows()-rowCount){
            cout << "[";
            for(int j=0; j<cols(); j++){
                if(j < colCount || j >= cols()-colCount)
                    cout << this->operator()(i,j) << "     ";
                else if(j == colCount)
                    cout << " ... ";
            }
            cout << "]" << endl;
        }else if(i == rowCount){
            cout << "..." << endl;
            cout << "..." << endl;
        }
    }
    cout << "]" << endl;
}

template<> void TensorD<3>::print3D(int depthCount, int rowCount, int colCount){
    printSize();
    std::cout << std::setprecision(3) << std::fixed;
    cout << "[";
    for(int d=0; d<depth(); d++){
        if(d < depthCount || d >= depth()-depthCount){
            printAtDepth(d, rowCount, colCount);
        }else if(d==depthCount){
            cout << "..." << endl;
            cout << "..." << endl;
        }
    }
    cout << "]" << endl;
}

template<> void TensorD<2>::print(){
    print2D();
}

template<> void TensorD<3>::print(){
    print3D();
}

std::string getSpaces(int spaces, bool empt = false){
    std::string s = "";
    if(empt) return s;
    for(int i=0; i<spaces; i++)
        s+=" ";
    return s;
}

void printEigen(const Eigen::MatrixXf& m, int colCount = 4, int rowCount = 4, int spaces = 4){
    std::cout << std::setprecision(3) << std::fixed;
    cout << "R: " << m.rows() << " C: " << m.cols() << endl;
    cout << "[";
    for(int i=0; i<m.rows(); i++){
        if(i < rowCount){
            cout << "[";
            for(int j=0; j<m.cols(); j++){
                bool empt = false;
                if(j == m.cols()-1) empt = true;
                if(j < colCount)
                    cout << m(i,j) << getSpaces(spaces,empt);
                else if(j == colCount)
                    cout << " ... ";
                else if(j >= m.cols()-colCount)
                    cout << m(i,j) << getSpaces(spaces,empt);
            }
            cout << "]" << endl;
        }else if(i == rowCount){
            cout << "..." << endl;
            cout << "..." << endl;
        }else if(i >= m.rows()-rowCount){
            cout << "[";
            for(int j=0; j<m.cols(); j++){
                bool empt = false;
                if(j == m.cols()-1) empt = true;
                if(j < colCount)
                    cout << m(i,j) << getSpaces(spaces,empt);
                else if(j == colCount)
                    cout << " ... ";
                else if(j >= m.cols()-colCount)
                    cout << m(i,j) << getSpaces(spaces,empt);
            }
            cout << "]" << endl;
        }
    }
    cout << "]" << endl;
}

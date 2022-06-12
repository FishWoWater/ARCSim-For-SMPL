#ifndef SMPL_HPP
#define SMPL_HPP

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigen>
// #include <jsoncpp/json/json.h>
// #include <jsoncpp/json/value.h>
// #include <jsoncpp/json/reader.h>
#include <json/json.h>
#include <json/value.h>
#include <json/reader.h>
#include "tensor.h"
#include <opencv2/opencv.hpp>
#include <fstream>

class SMPL{
public:
    Eigen::MatrixXf mV, mVTemp1, mVTemp2;
    Eigen::MatrixXf mF;
    Eigen::MatrixXf mPose;
    Eigen::MatrixXf mKintreeTable;
    Eigen::MatrixXf mJ, mJTemp1, mJTemp2;
    Eigen::MatrixXf mTrans;
    Eigen::MatrixXf mWeights;
    Eigen::MatrixXf mWeightsT;
    Eigen::MatrixXf vertSymIdxs;
    Eigen::MatrixXf mBetas;
    Eigen::MatrixXf mPose_flat;
    Eigen::MatrixXf mJR;
    typedef Eigen::Matrix<float, 4, 24> BlockMatrix;
    std::vector<Eigen::MatrixXf> weightedBlockMatrix1;
    std::vector<BlockMatrix> blocks;
    TensorD<3> mShapedDirsTensor;
    TensorD<3> mPoseDirsTensor;

    SMPL(){
        blocks.resize(4);
        weightedBlockMatrix1.resize(4);
        mBetas.resize(10,1);
        mPose.resize(24, 3); 
        mBetas.setZero();
        mPose_flat.resize(207, 1);
    }

    enum Part
    {
        BODY,   // 0
        LLEG,   // 1
        RLEG,   // 2
        LTORSO,    // 3
        LKNEE,  // 4
        RKNEE,  // 5
        MTORSO,  // 6
        LFOOT,  // 7
        RFOOT,  // 8
        UTORSO, // 9
        LLFOOT, // 10
        RRFOOT, // 11
        NECK,   // 12
        LSHOULDER,  // 13
        RSHOULDER,  // 14
        HEAD,   // 15
        LSHOULDER2, // 16
        RSHOULDER2, // 17
        LELBOW, // 18
        RELBOW, // 19
        LWRIST, // 20
        RWRIST, // 21
        LFINGERS, // 22
        RFINGERS, // 23
        TRANS, // 24
    };

    enum Shape
    {
        S0,
        S1,
        S2,
        S3,
        S4,
        S5,
        S6,
        S7,
        S8,
        S9
    };

    void setAllPoses(Eigen::VectorXf vec); 

    void setAllShapes(Eigen::VectorXf vec); 

    bool loadTensorFromJSON(const Json::Value& json, TensorD<3>& t, bool debug = false); 

    bool loadSparseFromJSON(const Json::Value& jsonval, Eigen::SparseMatrix<float>& m, int rows, int cols, bool debug = false); 

    bool loadPoseFromJSONFile(std::string filePath);

    bool loadModelFromJSONFile(std::string filePath);

    Eigen::Matrix4f rod(const Eigen::VectorXf& v, const Eigen::VectorXf& t); 

    bool updateModel(bool jointsOnly = false);

    bool loadEigenVecFromJSON(const Json::Value &json, vector<Eigen::MatrixXf> &t, bool debug = false);

    bool loadEigenFromJSON(const Json::Value &jsonval, Eigen::MatrixXf &m, bool neg = false);

    void saveToOBJ(const std::string);
};

#endif

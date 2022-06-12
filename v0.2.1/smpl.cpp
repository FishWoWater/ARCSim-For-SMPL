#include "smpl.h"

void SMPL::setAllPoses(Eigen::VectorXf vec){
        int row, col; 
        for(int i=0;i<72;i++){
            row = int(i / 3); 
            col = i % 3; 
            mPose(row, col) = vec(i);
        }
    }

void SMPL::setAllShapes(Eigen::VectorXf vec){
        for(int i=0;i<10;i++){
            mBetas(i, 0) = vec(i); 
        }
    }

void SMPL::saveToOBJ(const std::string savepath){
    std::ofstream fout;
    fout.open(savepath);
    for(int i=0;i<mVTemp2.rows();i++){
        fout << "v " << mVTemp2(i, 0) << " " << mVTemp2(i, 1) << " " << mVTemp2(i, 2) << std::endl;
    }
    std::cout << "saving " << mF.rows() << " faces" << std::endl;
    for(int i=0;i<mF.rows();i++){
        fout << "f " << mF(i, 0) + 1 << " " << mF(i, 1) + 1 << " " << mF(i, 2) + 1 << std::endl;
    }
    std::cout << "save mesh to" << savepath << std::endl;
}

bool SMPL::loadTensorFromJSON(const Json::Value& json, TensorD<3>& t, bool debug){
    int depth = json.size();
    int rows = json[0].size();
    int cols = json[0][0].size();
    if(debug){
        cout << "D: " << depth;
        cout << " R: " << rows;
        cout << " C: " << cols << endl;
    }

    t.resize({depth,rows,cols});

    for(int d=0; d<depth; d++){
        for(int r=0; r<rows; r++){
            for(int c=0; c<cols; c++){
                t(d,r,c) = json[d][r][c].asFloat();
            }
        }
    }

    return true;
}

bool SMPL::loadEigenVecFromJSON(const Json::Value& json, std::vector<Eigen::MatrixXf>& t, bool debug){
    int depth = json.size();
    int rows = json[0].size();
    int cols = json[0][0].size();
    if(debug){
        cout << "D: " << depth;
        cout << " R: " << rows;
        cout << " C: " << cols << endl;
    }

        t.resize(depth);
        for(Eigen::MatrixXf& m : t){
            m.resize(rows, cols);
        }

        for(int d=0; d<depth; d++){
            Eigen::MatrixXf& m = t[d];
            for(int r=0; r<rows; r++){
                for(int c=0; c<cols; c++){
                    t[d](r,c) = json[d][r][c].asFloat();
                }
            }
        }

        return true;
    }

bool SMPL::loadEigenFromJSON(const Json::Value& jsonval, Eigen::MatrixXf& m, bool neg){
        // Set Shape
        int rows = jsonval.size();
        if(!rows) { cerr << "Matrix Has no Rows" << endl; return false;}
        int cols = jsonval[0].size();
        if(cols == 0){
            m.resize(rows, 1);
            for(int i=0; i<rows; i++){
                m(i,0) = jsonval[i].asFloat();
            }
            return true;
        }
        if(rows == 0) rows = 1;
        m.resize(rows, cols);

        // Load Data
        if(rows > 1){
            for(int i=0; i<rows; i++){
                for(int j=0; j<cols; j++){
                    m(i,j) = jsonval[i][j].asFloat();
                    if(neg) m(i,j)*=-1;
                }
            }
        }else{
            throw std::runtime_error("Something wrong");
        }

        return true;
    }

bool SMPL::loadSparseFromJSON(const Json::Value& jsonval, Eigen::SparseMatrix<float>& m, int rows, int cols, bool debug){
        Eigen::MatrixXf intern;
        loadEigenFromJSON(jsonval, intern);
        std::vector<Eigen::Triplet<float>> tripletList;
        tripletList.reserve(intern.rows());
        for(int r=1; r<intern.rows(); r++){
            tripletList.push_back(Eigen::Triplet<float>(intern(r,0),intern(r,1),intern(r,2)));
            //cout << intern(r,1) << endl;
        }
        m = Eigen::SparseMatrix<float>(rows,cols);
        m.setFromTriplets(tripletList.begin(), tripletList.end());
        return true;
    }

bool SMPL::loadPoseFromJSONFile(std::string filePath){
    ifstream in(filePath);
    Json::Value root;
    in >> root;
    if(!root.size()){
        cerr << "Failed to load pose file" << endl;
        return false;
    }

    loadEigenFromJSON(root["pose"], mPose);
    loadEigenFromJSON(root["betas"], mBetas);
    loadEigenFromJSON(root["trans"], mTrans);
}

bool SMPL::loadModelFromJSONFile(std::string filePath){
    ifstream in(filePath);
    Json::Value root;
    in >> root;
    if(!root.size()){
        cerr << "Failed to load model file" << endl;
        return false;
    }

    loadEigenFromJSON(root["pose"], mPose);
    loadTensorFromJSON(root["shapedirs"], mShapedDirsTensor);
    loadTensorFromJSON(root["posedirs"], mPoseDirsTensor);
    loadEigenFromJSON(root["f"], mF);
    loadEigenFromJSON(root["kintree_table"], mKintreeTable);
    loadEigenFromJSON(root["J"], mJ);
    mJTemp1 = mJ;
    mJTemp2 = mJ;

    loadEigenFromJSON(root["trans"], mTrans);

    loadEigenFromJSON(root["v_posed"], mV);
    mV.conservativeResize(mV.rows(),mV.cols()+1);
    mV.col(mV.cols()-1) = Eigen::VectorXf::Ones(mV.rows());
    mVTemp1 = mV;
    mVTemp2 = mV;
    for(Eigen::MatrixXf& w : weightedBlockMatrix1) w.resize(4,mV.rows());

    loadEigenFromJSON(root["weights"], mWeights);
    mWeightsT = mWeights.transpose();

    loadEigenFromJSON(root["vert_sym_idxs"], vertSymIdxs);
    loadEigenFromJSON(root["J_regressor"], mJR);

    return true;
    }

Eigen::Matrix4f SMPL::rod(const Eigen::VectorXf& v, const Eigen::VectorXf& t){
        Eigen::Matrix4f m;
        cv::Mat src(cv::Size(1,3),CV_32FC1,cv::Scalar(0));
        src.at<float>(0) = v(0);
        src.at<float>(1) = v(1);
        src.at<float>(2) = v(2);
        cv::Mat dst;
        cv::Rodrigues(src, dst);
        m(0,0) = dst.at<float>(0,0);
        m(0,1) = dst.at<float>(0,1);
        m(0,2) = dst.at<float>(0,2);
        m(0,3) = t(0);
        m(1,0) = dst.at<float>(1,0);
        m(1,1) = dst.at<float>(1,1);
        m(1,2) = dst.at<float>(1,2);
        m(1,3) = t(1);
        m(2,0) = dst.at<float>(2,0);
        m(2,1) = dst.at<float>(2,1);
        m(2,2) = dst.at<float>(2,2);
        m(2,3) = t(2);
        m(3,0) = 0;
        m(3,1) = 0;
        m(3,2) = 0;
        m(3,3) = 1;
        return m;
    }

bool SMPL::updateModel(bool jointsOnly){

        // Create parent link table
        // {1: 0, 2: 0, 3: 0, 4: 1, 5: 2, 6: 3, 7: 4, 8: 5, 9: 6, 10: 7, 11: 8, 12: 9, 13: 9, 14: 9, 15: 12, 16: 13, 17: 14, 18: 16, 19: 17, 20: 18, 21: 19, 22: 20, 23: 21}
        std::map<int,int> parent;
        for(int i=1; i<mKintreeTable.cols(); i++){
            int key = mKintreeTable(1,i); int val = mKintreeTable(0,i);
            parent[key] = val;

//            cout << "{" << key << "," << val << "},";
        }

        std::vector<Eigen::Matrix4f> globalTransforms(24);
        std::vector<Eigen::Matrix4f> transforms(24);

        // Shape
        TensorD<3> AB = mShapedDirsTensor.dot(mBetas);
        for(int i=0; i<mVTemp1.rows(); i++){
            mVTemp1(i,0) = mV(i,0) + AB(i,0,0);
            mVTemp1(i,1) = mV(i,1) + AB(i,1,0);
            mVTemp1(i,2) = mV(i,2) + AB(i,2,0);
            mVTemp1(i,3) = 1;
        }
        // obtain poses_flat(207-dim vector from mPoses)
        for(int i=0; i<23; i++){
            cv::Mat src(cv::Size(1,3),CV_32FC1,cv::Scalar(0));
            src.at<float>(0) = mPose(i+1, 0);
            src.at<float>(1) = mPose(i+1, 1);
            src.at<float>(2) = mPose(i+1, 2);
            cv::Mat dst;
            cv::Rodrigues(src, dst);
            for(int j=0; j<3; j++){
                for(int k=0; k<3; k++){
                    if(j == k)  mPose_flat(i * 9 + 3 * j + k, 0) = dst.at<float>(j, k) - 1;
                    else mPose_flat(i * 9 + 3 * j + k, 0) = dst.at<float>(j, k);
                }
            }
        }


        // Shape J
        mJTemp2.row(0) = mJ.row(0);
        mJTemp1 = mJR * mVTemp1;

        TensorD<3> pose_blend_shapes = mPoseDirsTensor.dot(mPose_flat);
        for(int i=0; i<mVTemp1.rows(); i++){
            mVTemp1(i,0) = mVTemp1(i,0) + pose_blend_shapes(i,0,0);
            mVTemp1(i,1) = mVTemp1(i,1) + pose_blend_shapes(i,1,0);
            mVTemp1(i,2) = mVTemp1(i,2) + pose_blend_shapes(i,2,0);
            mVTemp1(i,3) = 1;
        }

        // cout << mJTemp1 << endl;

        // Body pose
        Eigen::Matrix4f& bodyPose = globalTransforms[0];
        bodyPose = rod(mPose.row(0), mJTemp1.row(0));

        // Global Transforms
        for(int i=1; i<globalTransforms.size(); i++){
            Eigen::Matrix4f& pose = globalTransforms[i];
            pose = globalTransforms[parent[i]] * rod(mPose.row(i), mJTemp1.row(i) - mJTemp1.row(parent[i]));
            mJTemp2(i,0) = pose(0,3);
            mJTemp2(i,1) = pose(1,3);
            mJTemp2(i,2) = pose(2,3);
        }

        // Trans
        for(int i=0; i<mJTemp2.rows(); i++){
            mJTemp2(i,0) = mJTemp2(i,0) + mTrans(0,0);
            mJTemp2(i,1) = mJTemp2(i,1) + mTrans(1,0);
            mJTemp2(i,2) = mJTemp2(i,2) + mTrans(2,0);
        }
        if(jointsOnly) return true;

        // Transforms
        for(int i=0; i<transforms.size(); i++){
            Eigen::Matrix4f& pose = transforms[i];
            Eigen::Vector4f jZero;
            jZero << mJTemp1(i,0), mJTemp1(i,1), mJTemp1(i,2), 0;
            Eigen::Vector4f fx = globalTransforms[i] * jZero; // Only apply rot to jVector
            Eigen::Matrix4f pack = Eigen::Matrix4f::Zero();
            pack(0,3) = fx(0);
            pack(1,3) = fx(1);
            pack(2,3) = fx(2);
            pose = globalTransforms[i] - pack; // Only minus t component from transform with rotated jVector
        }

        // Generate transform from weights
        for(int b=0; b<4; b++){
            BlockMatrix& block = blocks[b];
            for(int i=0; i<24; i++){
                block.col(i) = transforms[i].row(b);
            }
            weightedBlockMatrix1[b] = block*mWeightsT; // Column x VSize ~2ms
        }

        // Transform vertices with weight matrix
        for(int b=0; b<4; b++){
            Eigen::MatrixXf& block = weightedBlockMatrix1[b];
            for(int i=0; i<mV.rows(); i++){
                mVTemp2(i,b) = mVTemp1.row(i) * block.col(i);
            }
        }

        // Final transform
        for(int i=0; i<mVTemp2.rows(); i++){
            mVTemp2(i,0) = mVTemp2(i,0) + mTrans(0,0);
            mVTemp2(i,1) = mVTemp2(i,1) + mTrans(1,0);
            mVTemp2(i,2) = mVTemp2(i,2) + mTrans(2,0);
        }

        return true;

}

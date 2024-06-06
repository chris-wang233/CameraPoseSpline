#include "SE3spline.h"
#include "sophus/so3.hpp"
#include <map>

class SO3Spline{
    public:
     void compute_spline_para(Matrix3d R1,Eigen::Matrix3d W0,Eigen::Matrix3d W1_,Matrix3d R0){
        //这个r1是R0T*R1
        //Eigen::Vector3d _w1;
        Matrix3d r1 = SO3_::log(R1);
        Matrix3d index = SO3_::log(R0);
        Eigen::Matrix3d w1;
        Vector3d w0 = SO3_::anti_hat(W0);
        Vector3d w1_ = SO3_::anti_hat(W1_);
        //Matrix3d rj = SO3_::right_jacobian(r1);
        //Matrix3d rj_in = rj.inverse();
        w1_ =SO3_::right_jacobian_inverse(index)*w1_;

        //w1_ =rj_in*w1_;

        //W0 = w0;
       // W1_ = SO3_::hat(w1_);
        spline_c = W0;
        spline_a = -2*r1+W0+W1_;
        spline_b = 3*r1-2*W0-W1_;

        //spline_a = -2*r1+W0+w1_;
        //spline_b = 3*r1-2*W0-w1_;

        R0_ = R0;
        //spline_a = W0 - r1;
       // spline_b = 2*r1 - W0;
        
    }

    void compute_spline_para(Vector3d Lie_alg,Eigen::Matrix3d W0,Eigen::Matrix3d W1_,Matrix3d R0){
        //这个r1是R0T*R1
        //Eigen::Vector3d _w1;
        Matrix3d r1 = SO3_::hat(Lie_alg);
        Eigen::Matrix3d w1;
        Vector3d w0 = SO3_::anti_hat(W0);
        Vector3d w1_ = SO3_::anti_hat(W1_);
        Matrix3d index = SO3_::log(R0);
        //Matrix3d rj = SO3_::right_jacobian(r1);
        //Matrix3d rj_in = rj.inverse();
        w1_ =SO3_::right_jacobian_inverse(index)*w1_;
        
        //W1_ = SO3_::hat(w1_);
        //w1_ =rj_in*w1_;

        //W0 = w0;
       // W1_ = SO3_::hat(w1_);
        spline_c = W0;
        spline_a = -2*r1+W0+W1_;
        spline_b = 3*r1-2*W0-W1_;

        R0_ = R0;
        //spline_a = W0 - r1;
       // spline_b = 2*r1 - W0;
        
    }

    void compute_spline_para(map<double,Vector3d> Lie_alg){

    }

    Matrix3d compute_rotation(double time_persent,Matrix3d R0){
        Matrix3d spline_algebra = spline_a*time_persent*time_persent*time_persent+spline_b*time_persent*time_persent+
        spline_c*time_persent;
        Vector3d sp_al= SO3_::anti_hat(spline_algebra);
       // Matrix3d delta_rotation;
/*         if(R0_.trace())
        {   Sophus::SO3d Rs(R0_);
            Sophus::Vector3d Sop_R =Rs.log();
            Vector3d Sop_R_ = Sop_R;
            Matrix3d index = SO3_::hat(sp_al+Sop_R_);
            delta_rotation = SO3_::exp(index);
        }
        else{
            delta_rotation =R0_*SO3_::exp(spline_algebra);
        } */
        //Matrix3d index =  R_lie+spline_algebra;
        
        //Matrix3d spline_algebra = spline_a*time_persent*time_persent+spline_b*time_persent;
        Matrix3d Ja =  SO3_::right_jacobian(SO3_::log(R0_));
        Matrix3d index = Ja*spline_algebra;
        Matrix3d delta_rotation = SO3_::exp(index);
        Matrix3d delta_rotation_ = SO3_::exp(spline_algebra);
        return R0_*delta_rotation_;
    }

    Matrix3d SVDSpline(double time_persent,double start_time,double end_time,
                        Matrix3d start_R,Matrix3d end_R,Vector3d start_v,Vector3d end_v){
        Matrix3d M3,M2,M1,M0;
        Matrix3d w0 = SO3_::hat(start_v);
        Matrix3d w1 = SO3_::hat(end_v);
        Matrix3d dR0 = w0*start_R;
        Matrix3d dR1 = w1*end_R;
        start_time = 0;
        end_time = 1;

        double delta_t = end_time-start_time;
        //为什么不是R0^T*R1？
        Matrix3d delta_x = end_R - start_R;
        //为什么不是在李代数空间进行计算？？与上面的方程有同样的问题，就是该方法对旋转矩阵直接进行加减法
        Matrix3d delta_v_log = SO3_::log(dR1)-SO3_::log(dR0);
        Matrix3d delta_v = dR1-dR0;

        M3 = 6*((dR0+dR1)/(delta_t*delta_t))-(12*(delta_x/(delta_t*delta_t*delta_t)));
        M2 = (delta_v/delta_t)-M3*(end_time+start_time)*0.5;
        M1 = dR0-M3*0.5*start_time*start_time-M2*start_time;
        M0 = start_R - M3*start_time*start_time*start_time/6-M2*start_time*start_time/2-M1*start_time;

        Matrix3d M_t = M3*time_persent*time_persent*time_persent+M2*time_persent*time_persent+M1*time_persent+M0;
        M_t = M_t*0.4*Eigen::Matrix3d::Identity();
        //printf("start time is:%f\n",start_time);
        //printf("end time is:%f\n",end_time);
        //printf("MT is:\n%f,%f,%f\n%f,%f,%f\n%f,%f,%f\n\n",M_t(0,0),M_t(0,1),M_t(0,2),M_t(1,0),M_t(1,1),M_t(1,2),M_t(2,0),M_t(2,1),M_t(2,2));
        JacobiSVD<Matrix3d> svd(M_t,ComputeFullU|ComputeFullV);
        Matrix3d R_t = svd.matrixU()*svd.matrixV().transpose();
        return R_t;
    }

    void get_para(Vector3d &W1,Vector3d &W2,Vector3d &W3)
    {
       Vector3d W1_ = SO3_::anti_hat(spline_a);
       Vector3d W2_ = SO3_::anti_hat(spline_b);
       Vector3d W3_ = SO3_::anti_hat(spline_c);

       W1 = W1_;
       W2 = W2_;
       W3 = W3_;
    }

    void updatePara(Vector3d W1,Vector3d W2,Vector3d W3){
        spline_a = SO3_::hat(W1);//P3
        spline_b = SO3_::hat(W2);//P2
        spline_c = SO3_::hat(W3);//P1
    }

    private:
    Eigen::Matrix3d spline_a;
    Eigen::Matrix3d spline_b;
    Eigen::Matrix3d spline_c;
    Eigen::Matrix3d R0_;

};

class SO3Spline_Least_Square_muti{
    public:
        SO3Spline_Least_Square_muti(std::map<double,Eigen::Vector3d> time_pose_,Eigen::Vector3d start_pose_,
        Eigen::Vector3d end_pose_,Eigen::Vector3d start_velocity_,Eigen::Vector3d end_velocity_):
        time_pose(time_pose_),start_pose(start_pose_),end_pose(end_pose_),start_velocity(start_velocity_),
        end_velocity(end_velocity_),para_factor(4,4+time_pose.size()){};

        void compute_para(){
            std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();
/*             P0=para_factor(0,0)*start_pose+para_factor(0,1)*end_pose+para_factor(0,2)*start_velocity+para_factor(0,3)*end_velocity+para_factor(0,4)*it->second;
            P1=para_factor(1,0)*start_pose+para_factor(1,1)*end_pose+para_factor(1,2)*start_velocity+para_factor(1,3)*end_velocity+para_factor(1,4)*it->second;
            P2=para_factor(2,0)*start_pose+para_factor(2,1)*end_pose+para_factor(2,2)*start_velocity+para_factor(2,3)*end_velocity+para_factor(2,4)*it->second;
            P3=para_factor(3,0)*start_pose+para_factor(3,1)*end_pose+para_factor(3,2)*start_velocity+para_factor(3,3)*end_velocity+para_factor(3,4)*it->second; */

            P0 =0*para_factor(0,0)*start_pose+para_factor(0,1)*end_pose+0.01*para_factor(0,2)*start_velocity+0.01*para_factor(0,3)*end_velocity;
            P1 =0*para_factor(1,0)*start_pose+para_factor(1,1)*end_pose+0.01*para_factor(1,2)*start_velocity+0.01*para_factor(1,3)*end_velocity;
            P2 =0*para_factor(2,0)*start_pose+para_factor(2,1)*end_pose+0.01*para_factor(2,2)*start_velocity+0.01*para_factor(2,3)*end_velocity;
            P3 =0*para_factor(3,0)*start_pose+para_factor(3,1)*end_pose+0.01*para_factor(3,2)*start_velocity+0.01*para_factor(3,3)*end_velocity;

            int i=1;
            for(std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();it!=time_pose.end();it++){
                P0 += para_factor(0,3+i)*it->second;
                P1 += para_factor(1,3+i)*it->second;
                P2 += para_factor(2,3+i)*it->second;
                P3 += para_factor(3,3+i)*it->second;
                i++;
/*             printf("P1 is: %f,%f,%f\n",P1.x(),P1.y(),P1.z());
            printf("P2 is: %f,%f,%f\n",P2.x(),P2.y(),P2.z());
            printf("P3 is: %f,%f,%f\n\n",P3.x(),P3.y(),P0.z()); */
            }
            //printf("end loop\n\n");
        }

        //利用最小二乘法计算矩阵系数，即(A^T*A)^(-1)*A^T
        void compute_para_factor(){
            
            Eigen::MatrixXd A(4+time_pose.size(),4);

/*             A<<
                1,0,0,0,
                1,1,1,1,
                0,1,0,0,
                0,1,2,3; */

            A(0,0)=1;
            A(0,1)=A(0,2)=A(0,3)=0;

            A(1,0)=A(1,1)=A(1,2)=A(1,3)=1;

            A(2,0)=A(2,2)=A(2,3)=0;
            A(2,1)=0.01;

            A(3,0)=0;
            A(3,1)=0.01;
            A(3,2)=0.02;
            A(3,3)=0.03;
           // printf("hey hey damm\n\n");
            int i =1;
            for(std::map<double,Eigen::Vector3d>::iterator it=time_pose.begin();it!=time_pose.end();it++){
            
                A(3+i,0)=1;
                A(3+i,1)=it->first;
                A(3+i,2)=it->first*it->first;
                A(3+i,3)=it->first*it->first*it->first;
                i++;
/*             A<<
                1,0,0,0,
                1,1,1,1,
                0,1,0,0,
                0,1,2,3,
                1,it->first,it->first*it->first,it->first*it->first*it->first; */


            }
           Eigen::Matrix4d output_first = A.transpose()*A;
           para_factor = output_first.completeOrthogonalDecomposition().pseudoInverse()*A.transpose();
            compute_para();

        }

        //赋值
        void get_para(Eigen::Vector3d &P0_,Eigen::Vector3d &P1_,Eigen::Vector3d &P2_,Eigen::Vector3d &P3_){
            P0_ = P0;
            P1_ = P1;
            P2_ = P2;
            P3_ = P3;

        }


    private:
        //map 键值为time_persent，内容为三维pose信息
        std::map<double,Eigen::Vector3d> time_pose;
        //端点信息
        Eigen::Vector3d start_pose;
        Eigen::Vector3d end_pose;
        Eigen::Vector3d start_velocity;
        Eigen::Vector3d end_velocity;
        Eigen::MatrixXd para_factor;
        Eigen::Vector3d P0;
        Eigen::Vector3d P1;
        Eigen::Vector3d P2;
        Eigen::Vector3d P3;
};
/*
FILE amp/rrt.cpp

LICENSE

Copyright 2017 Samantha Stoneman
Research Associate, German Aerospace Center (DLR)

This file is part of the Articulated-robot Motion Planner (Amp) library.

Amp is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Amp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Amp. If not, see <http://www.gnu.org/licenses/>.
*/

#include "amp/rrt.hpp"

namespace amp {


template<typename scalar>
rrt<scalar>
::rrt(const std::shared_ptr<planner_t>& planner_in,
      scalar r_radius_max_in,
      scalar q_radius_max_in) :
    planner_(planner_in),
    kdtree_(kd_create(7)),
    search_type(DEFAULT),
    generator(std::chrono::high_resolution_clock::now().time_since_epoch().count()),
    radius_max(Eigen::Matrix<scalar,2,1>::Zero())
{
    seed_rng();
    radius_max(0) = r_radius_max_in;
    radius_max(1) = q_radius_max_in;
}


template<typename scalar>
rrt<scalar>
::~rrt()
{
    kd_free(kdtree_);
}


template<typename scalar>
void rrt<scalar>
::seed_rng()
{
    for(int didx=0; didx<1; ++didx) {
        if (didx!=2)
            distributions.push_back(uni_t(-1.*planner_->tcp_bound, planner_->tcp_bound));
        else
            distributions.push_back(uni_t(0., planner_->tcp_bound));
        distributions.push_back(uni_t(-1.,1.)); // quaternion limits
    }
    // randomly sample time at waypoint
    //uni_t time_distribution = uni_t(5., final_time_des-20.);
    //node_sample->time_curr = time_distribution(generator);
}


template<typename scalar>
typename rrt<scalar>::rrt_result rrt<scalar>
::run(const std::shared_ptr<const node_t>& root_, const int stop_iter, const bool connect_to_goal_first)
{
#ifdef _VERBOSE_
    std::cerr << "--Running RRT-- " << std::endl;
#endif //_VERBOSE_
    std::shared_ptr<branches_t> root_branch = create_branch(root_);
    insert_node(root_branch);

    int rrt_iter = 0;
    bool connected = false;
    bool sample_connected;
    if (connect_to_goal_first)
        sample_connected=true;
    else
        sample_connected=false;

    std::shared_ptr<traj_t> traj_;
    std::shared_ptr<const node_t> node_last;
    std::shared_ptr<node_t> node_sample;
    std::shared_ptr<const node_t> node_nearest;
    branches_t* branch_nearest;
    std::shared_ptr<branches_t> branch_last;
    typename branches_t::connections_t::const_iterator iter;

    while (rrt_iter < stop_iter && connected==false) {
        if (sample_connected) {
            node_last = tree_.back()->node_;
            node_sample = std::make_shared<node_t>(goal_);
            node_sample->time_curr = 2. + node_last->time_curr;
#ifdef _VERBOSE_
            std::cerr << "Node From: ";
            std::cerr << node_last->state.block(0,0,7,1).transpose() << std::endl;
            std::cerr << "Node To: ";
            std::cerr << node_sample->state.block(0,0,7,1).transpose() << std::endl;
#endif //_VERBOSE_
            if (planner_t::SUCCESS == connect(traj_,node_last,node_sample)) {
                traj_->node(node_sample,traj_->via_points()-1);
                branch_last = tree_.back();
                branch_last->insert_branch(traj_);
                iter = --branch_last->outgoing_branches_.end();
                std::shared_ptr<branches_t> branch_goal =
                        std::make_shared<branches_t>(node_sample,branch_last.get(),iter);
                insert_node(branch_goal);
                solutions.push_back(branch_goal);
                connected = true;
            }
#ifdef _LOG_ALL_DATA_
            else
                log_failed(node_sample, traj_);
#endif //_LOG_ALL_DATA_
        }
        if (!connected) {
            node_sample = create_sample();
            branch_nearest = nearest_node(node_sample);
            node_nearest = branch_nearest->node_;
            node_sample->time_curr = 2. + node_nearest->time_curr;
            shorten_sample_distance(node_sample,node_nearest);
#ifdef _VERBOSE_
            std::cerr << "\nIteration " << rrt_iter << ": Tree size... " << tree_.size();
            std::cerr << ", Make a new sample..." << std::endl;
            std::cerr << "Node From: ";
            std::cerr << node_nearest->state.block(0,0,7,1).transpose() << std::endl;
            std::cerr << "Node To: ";
            std::cerr << node_sample->state.block(0,0,7,1).transpose() << std::endl;
#endif //_VERBOSE_
            if(planner_t::SUCCESS == connect(traj_, node_nearest, node_sample)) {
                traj_->node(node_sample,traj_->via_points()-1);
                branch_nearest->insert_branch(traj_);
                iter = --branch_nearest->outgoing_branches_.end();
                std::shared_ptr<branches_t> branch_sample =
                        std::make_shared<branches_t>(node_sample,branch_nearest,iter);
                insert_node(branch_sample);
                sample_connected = true;
            }
            else {
                sample_connected = false;
#ifdef _LOG_ALL_DATA_
                log_failed(node_sample, traj_);
#endif //_LOG_ALL_DATA_
            }
        }
        ++rrt_iter;
    }

#ifdef _LOG_ALL_DATA_
    write_data_to_file(samples_failed,trajs_failed,"failed");
#endif //_LOG_ALL_DATA_
    rrt_result res;
    res.iterations = rrt_iter;
    res.edges = tree_.size();
    res.goal_connected = connected;
    return(res);
}


template<typename scalar>
void rrt<scalar>
::shorten_sample_distance(const std::shared_ptr<node_t>& sample_in,
                          const std::shared_ptr<const node_t>& nearest_in)
{
    Eigen::Matrix<scalar,2,1> sample_radius = Eigen::Matrix<scalar,2,1>::Zero();

    sample_radius(0) = (sample_in->state.block(0,0,3,1) -
                        nearest_in->state.block(0,0,3,1)).norm();
    sample_radius(1) = (sample_in->state.block(3,0,4,1) -
                        nearest_in->state.block(3,0,4,1)).norm();

    if((sample_radius.array()>radius_max.array()).any()) {
        sample_in->state.block(0,0,3,1) =
                nearest_in->state.block(0,0,3,1) +
                ((sample_in->state.block(0,0,3,1) -
                  nearest_in->state.block(0,0,3,1)).normalized() * radius_max(0));
        sample_in->state.block(3,0,4,1) =
                nearest_in->state.block(3,0,4,1) +
                ((sample_in->state.block(3,0,4,1) -
                  nearest_in->state.block(3,0,4,1)).normalized() * radius_max(1));
    }
}


template<typename scalar>
void rrt<scalar>
::log_failed(const std::shared_ptr<const node_t>& sample_failed,
             const std::shared_ptr<const traj_t>& traj_failed)
{
    samples_failed.push_back(sample_failed);
    trajs_failed.push_back(traj_failed);
}



template<typename scalar>
void rrt<scalar>
::write_data_to_file(sample_list_t& samples, traj_list_t& trajs, const char* filename)
{
    int idx=0;
    matrix_t sample_nodes = matrix_t::Zero(samples.size(),8);
    for(auto sample=samples.begin(); sample!=samples.end(); ++sample, ++idx) {
        std::shared_ptr<const node_t> sample_ = *sample;
        sample_nodes(idx,0) = sample_->time_curr;
        sample_nodes.block(idx,1,1,7) = sample_->state.block(0,0,7,1).transpose();
    }
    rwm::write_matrix_to_file(rwm::generic_name(filename,"samples",0).c_str(),sample_nodes);

    idx = 0;
    for(auto traj=trajs.begin(); traj!=trajs.end(); ++traj, ++idx) {
        std::shared_ptr<const traj_t> traj_ = *traj;
        write_traj_to_file(traj_,filename,idx);
    }
}


template<typename scalar>
void rrt<scalar>
::write_traj_to_file(const std::shared_ptr<const traj_t>& traj_, const char* filename, const int idx)
{
    matrix_t state_and_time = matrix_t::Zero(traj_->via_points(), planner_->dof()+1);
    state_and_time.block(0,0,state_and_time.rows(),1) = traj_->times();
    state_and_time.block(0,1,state_and_time.rows(),planner_->dof()) = traj_->joint_dimensions();
    rwm::write_matrix_to_file(rwm::generic_name(filename,"states",idx).c_str(),state_and_time);

    matrix_t dstate_and_time = matrix_t::Zero(traj_->via_points(), planner_->dof()+1);
    dstate_and_time.block(0,0,dstate_and_time.rows(),1) = traj_->times();
    dstate_and_time.block(0,1,dstate_and_time.rows(),planner_->dof()) = traj_->joint_dimensions(1);
    rwm::write_matrix_to_file(rwm::generic_name(filename,"dstates",idx).c_str(),dstate_and_time);

    matrix_t tcp_and_time = matrix_t::Zero(traj_->via_points(), 8);
    tcp_and_time.block(0,0,tcp_and_time.rows(),1) = traj_->times();
    tcp_and_time.block(0,1,tcp_and_time.rows(),7) = traj_->dimensions(0);
    rwm::write_matrix_to_file(rwm::generic_name(filename,"tcps",idx).c_str(),tcp_and_time);

    matrix_t dtcp_and_time = matrix_t::Zero(traj_->via_points(), 8);
    dtcp_and_time.block(0,0,dtcp_and_time.rows(),1) = traj_->times();
    dtcp_and_time.block(0,1,dtcp_and_time.rows(),7) = traj_->dimensions(1);
    rwm::write_matrix_to_file(rwm::generic_name(filename,"dtcps",idx).c_str(),dtcp_and_time);

    matrix_t manip_and_time = matrix_t::Zero(traj_->via_points(), 2);
    manip_and_time.block(0,0,manip_and_time.rows(),1) = traj_->times();
    manip_and_time.block(0,1,manip_and_time.rows(),1) = traj_->manipulability();
    rwm::write_matrix_to_file(rwm::generic_name(filename,"manips",idx).c_str(),manip_and_time);
}


template<typename scalar>
void rrt<scalar>
::parse_solutions()
{
    int idx = 0;
    for(auto goal_curr=solutions.begin(); goal_curr!=solutions.end(); ++goal_curr, ++idx) {
        const branches_t* branch_curr = (*goal_curr).get();
        std::shared_ptr<const node_t> node_curr = branch_curr->node_;
        const std::shared_ptr<const node_t> node_initial = tree_.front()->node_;

        int n_steps_rrt_trajectory = 1;
        matrix_t rrt_trajectory;

        std::shared_ptr<const traj_t> traj_prev;
        std::list<std::shared_ptr<const traj_t> > intermediate_trajs;

        while (node_curr!=node_initial) {
            const branches_t* branch_prev = branch_curr->incoming_branches.first;
            traj_prev = *(branch_curr->incoming_branches.second);
            intermediate_trajs.insert(intermediate_trajs.begin(),traj_prev);
            node_curr = branch_prev->node_;
            branch_curr = branch_prev;
            n_steps_rrt_trajectory += intermediate_trajs.front()->via_points()-1 ;
        }

        if (intermediate_trajs.size() > 1) {
            int count_traj = 1;
            int count_via = 0;

            rrt_trajectory =
                    matrix_t::Zero(n_steps_rrt_trajectory, intermediate_trajs.front()->n_columns());
            auto traj_curr = intermediate_trajs.begin();
            count_via = (*traj_curr)->via_points();

            rrt_trajectory.block(0,0,(*traj_curr)->via_points(),(*traj_curr)->n_columns()) =
                    (*traj_curr)->whole_trajectory();
            ++traj_curr;
            for (; traj_curr!=intermediate_trajs.end(); ++traj_curr) {
                rrt_trajectory.block(count_via, 0, (*traj_curr)->via_points()-1, (*traj_curr)->n_columns()) =
                        (*traj_curr)->whole_trajectory().block(1, 0, (*traj_curr)->via_points()-1, (*traj_curr)->n_columns());
                count_via += (*traj_curr)->via_points()-1;
                ++count_traj;
            }
            std::shared_ptr<traj_t> solution_curr =
                    std::make_shared<traj_t>(intermediate_trajs.back()->dofj,
                                             intermediate_trajs.back()->ders,
                                             rrt_trajectory.rows());
            solution_curr->set_data(rrt_trajectory);
            trajs_solution.push_back(solution_curr);
        }
    }
}


template<typename scalar>
void rrt<scalar>
::parse_tree()
{
    /** start at last branch and iterate through tree */
    int idx = 0;
    idx = tree_.size()-1;
    auto branch=--tree_.end();
    for(; branch!=tree_.begin(); --idx) {
        std::shared_ptr<const branches_t> branch_curr = (*branch);
        trajs_succeeded.push_back(*(branch_curr->incoming_branches.second));
        samples_succeeded.push_back(branch_curr->node_);
        --branch;
    }
    samples_succeeded.push_back((*branch)->node_);
    write_data_to_file(samples_succeeded,trajs_succeeded, "tree");
}


template<typename scalar>
std::shared_ptr<typename rrt<scalar>::node_t>
rrt<scalar>::create_sample()
{
    std::shared_ptr<node_t> sample = std::make_shared<node_t>(planner_->dof(),3);
    //positions
    for (int idx=0; idx<3; ++idx) {
        sample->state(idx,0) = distributions[0](generator);
    }
    //quaternions
    for (int idx=3; idx<7; ++idx) {
        sample->state(idx,0) = distributions[1](generator);
    }
    return(sample);
}


template<typename scalar>
std::shared_ptr<typename rrt<scalar>::branches_t>
rrt<scalar>::create_branch(const std::shared_ptr<const node_t>& root_in)
{
    std::shared_ptr<branches_t> root_new = std::make_shared<branches_t>(root_in);
    return(root_new);
}


template<typename scalar>
void rrt<scalar>::insert_node(const std::shared_ptr<branches_t>& branch_in)
{
    tree_.push_back(branch_in);
    const scalar* position_ptr =  branch_in->node_->state.data(); // colmajor dof x ders matrix!
    void* branch_ptr = reinterpret_cast<void*>(branch_in.get());
    kd_insert(kdtree_,position_ptr,branch_ptr);
}


template<typename scalar>
typename rrt<scalar>::branches_t*
rrt<scalar>::nearest_node(const std::shared_ptr<node_t>& node_in)
{
    const scalar* position_ptr =  node_in->state.data(); // colmajor dof x ders matrix!
    void* branch_ptr = kd_res_item_data(kd_nearest(kdtree_,position_ptr));
    branches_t* branch_raw_ptr = reinterpret_cast<branches_t*>(branch_ptr);
    return(branch_raw_ptr);
}


/* template<typename scalar>
typename rrt<scalar>::rrt::tree_t
rrt<scalar>::nearest_nodes_in_radius(std::shared_ptr<node_t> node_in, const scalar)
{

}
*/


template<typename scalar>
typename rrt<scalar>::planner_t::return_type_t
rrt<scalar>::connect(std::shared_ptr<traj_t>& traj_result,
                     const std::shared_ptr<const node_t>& node_from,
                     const std::shared_ptr<const node_t>& node_to)
{
    // 5 is the # of free params
    planner_->set_desired_tcp_trajectory(5,node_from,node_to);
    typename planner_t::return_type_t result =
            planner_->run_adaptive(node_from->joints.block(0,0,planner_->dof(),1),
                                   node_from->joints.block(0,1,planner_->dof(),1));
    traj_result = planner_->tcp_traj();
    return(result);
}


// instantiate robot class
template class rrt<double>;

} //end namespace amp

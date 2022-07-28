classdef MH_genetic_it_data
    properties
        %% ints_
        leader_count;
        gen_count;
        nsuccess;
        ncandidates;
        bleader_rid;
        wleader_rid;
        nreplace;
        ndup;
        ndup_unique;
        nredraw;
        nreload;
        Class_count;
        event_block_count;
        %% dubs_
        rho2;
        br2;
        wr2;
        prob_best;
        prob_worst;
    end
    methods
        function obj = MH_genetic_it_data(ints_,dubs_)
            obj.leader_count=ints_(1);
            obj.gen_count=ints_(2);
            obj.nsuccess=ints_(3);
            obj.ncandidates=ints_(4);
            obj.bleader_rid=ints_(5);
            obj.wleader_rid=ints_(6);
            obj.nreplace=ints_(7);
            obj.ndup=ints_(8);
            obj.ndup_unique=ints_(9);
            obj.nredraw=ints_(10);
            obj.nreload=ints_(11);
            obj.Class_count=ints_(12);
            obj.event_block_count=ints_(13);

            obj.rho2=dubs_(1);
            obj.br2=dubs_(2);
            obj.wr2=dubs_(3);
            obj.prob_best=dubs_(4);
            obj.prob_worst=dubs_(5);
        end
    end
end

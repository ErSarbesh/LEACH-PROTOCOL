%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               LEACH                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
warning off
xm=100;      %diameters of sensor network
ym=100;
sink.x=100;  %distance of base station from the network
sink.y=100;
n = 100;         %no of nodes
p=0.1;          %probibilty of a node to become cluster head
Eo=0.5;          %energy supplied to each node
ETX=50*0.000000001;     %transmiter energy per node
ERX=50*0.000000001;        %reciever energy per mode
Efs=10*0.0000000000001;     %amplification energy when d is less than d0
Emp=0.0013*0.000000000001;      %amplification energy  when d is greater than d0
%Data Aggregation Energy
EDA=5*0.000000001;
rmax=3500;           %no of rounds
do=sqrt(Efs/Emp);       %distance between cluster head and base station

for i=1:1:n
    S(i).xd=rand(1,1)*xm;         %it will distribute the nodes in 1 dimension in x axis randomly.
    S(i).yd=rand(1,1)*ym;           %it will distribute the nodes in 1 dimension in y axis randoml.
    S(i).G=0;                        % as the no of node that have been cluster head is zero 0
    S(i).E=Eo%%*(1+rand*a);                %?
    %initially there are no cluster heads only nodes
    S(i).type='N';
end
S(n+1).xd=sink.x;   %assume that base station is also a node sp total no of nodes is n and with base station  it is n+1
S(n+1).yd=sink.y;
countCHs=0;         %the number of Stateflow objects in the current context.
cluster=1;              %first cluster is selected
flag_first_dead=0;
flag_teenth_dead=0;
flag_all_dead=0;
dead=0;
first_dead=0;
teenth_dead=0;
all_dead=0;
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
% counter for sleep nodes
s=0;

for r=0:1:rmax
    r
    if(mod(r, round(1/p) )==0) %remainder
        for i=1:1:n
            S(i).G=0;            % it will assign to the nodes that have not been cluster head .
        end
    end
    dead=0;
    for i=1:1:n
        
        if (S(i).E<=0)
            dead=dead+1;
            
            if (dead==1)
                if(flag_first_dead==0)
                    first_dead=r;
                    flag_first_dead=1;
                end
            end
            
            if(dead==0.1*n)
                if(flag_teenth_dead==0)
                    teenth_dead=r;
                    flag_teenth_dead=1;
                end
            end
            if(dead==n)
                if(flag_all_dead==0)
                    all_dead=r;
                    flag_all_dead=1;
                end
            end
        end
        if (S(i).E>0)
            S(i).type='N';
        end
    end
    STATISTICS.DEAD(r+1)=dead;
    STATISTICS.ALLIVE(r+1)=allive-dead;
    
    countCHs=0;
    cluster=1;
%     CH selection
    for i=1:1:n
        if(S(i).E>=0)  %Checking threshold and nmbr of sleep nodes
            temp_rand=rand;
            if ( (S(i).G)<=0)
                
                if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                    countCHs=countCHs+1;
                    packets_TO_BS=packets_TO_BS+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                    S(i).type='C';
                    S(i).G=round(1/p)-1;
                    C(cluster).xd=S(i).xd;
                    C(cluster).yd=S(i).yd;
                    distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                    C(cluster).distance=distance;
                    C(cluster).id=i;
                    X(cluster)=S(i).xd;
                    Y(cluster)=S(i).yd;
                    cluster=cluster+1;
                    distance;
                    if (distance>do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*(distance*distance*distance*distance ));
                    end
                    if (distance<=do)
                        S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*(distance * distance ));
                    end
                end 
            end
        end
    end
    STATISTICS.COUNTCHS(r+1)=countCHs;
    %   Association of nodes
    for i=1:1:n
        if ( S(i).type=='N' && S(i).E>0)
            if(cluster-1>=1)
                min_dis=Inf;
                min_dis_cluster=0;
                for c=1:1:cluster-1
                    temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster=c;
                    end
                end
                
                
                min_dis;
                if (min_dis>do)
                    S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis *min_dis * min_dis * min_dis));
                end
                if (min_dis<=do)
                    S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
                end
                    S(C(min_dis_cluster).id).E =S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 );
                    packets_TO_CH=packets_TO_CH+1;
                
                    S(i).min_dis=min_dis;
                    S(i).min_dis_cluster=min_dis_cluster; 
     else
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
           if (min_dis>do)
               S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis *min_dis * min_dis * min_dis));
           end
           if (min_dis<=do)
               S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis));
           end
           packets_TO_BS=packets_TO_BS+1;
            end
        end
    end    
                
    STATISTICS.PACKETS_TO_CH(r+1)=packets_TO_CH;
    STATISTICS.PACKETS_TO_BS(r+1)=packets_TO_BS;
    

end

 
figure(1);
r=0:rmax;
plot(r,STATISTICS.ALLIVE(r+1),'-b')%, r,STATISTICS.ALLIVE2(r+1),'-g', r,STATISTICS.ALLIVE3(r+1),'-y')
legend('LEACH');%'Z-SEP  m=0.1,a=2');
xlabel('No. of Rounds (r)');
ylabel('No. of  Nodes Allive');

figure(2);
plot(r, STATISTICS.DEAD(r+1),'r')
legend('LEACH');
xlabel('Rounds');
ylabel('Number of Dead Nodes');
title('Number of Dead Nodes vs Number of Rounds');

% figure(3)
% for r=1:rmax+1
%     En(r)=0;
%     for i=1:n
%         if (S(i).E>0)
%             En(r)=En(r)+S(i).E;
%         end
%     end
%     avg_residual_energy(r)=En(r)/n;
% end
% plot(avg_residual_energy)
% legend('LEACH');
% xlabel('Rounds');
% ylabel('Average Residual Energy');
% title('Average Residual Energy vs Number of Rounds');
% 
% figure(4)
% plot(r,PACKETS_TO_BS,'-b')%r,STATISTICS.DEAD1(r+1),'-r');%,r,STATISTICS.DEAD2(r+1),'-g')
% legend('LEACH');%'Z-SEP  m=0.1,a=2');
% xlabel('Rounds');
% ylabel('Throughput');


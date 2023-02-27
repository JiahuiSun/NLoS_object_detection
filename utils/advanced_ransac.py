import scipy.io as io
import math
import matplotlib
import matplotlib.pyplot as plt
import random
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cmx
import copy

font2 = {'family': 'Times New Roman',
        'weight': 'bold',
         'size': 16,
         }

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.weight'] = 'bold'

def kalman(z_measure,x_last=0,p_last=0,Q=0.018,R=0.0542):
    x_mid = x_last
    p_mid = p_last + Q
    kg = p_mid/(p_mid + R)
    x_now = x_mid + kg*(z_measure - x_mid)
    p_now = (1-kg)*p_mid
    p_last = p_now
    x_last = x_now
    return x_now,p_last,x_last

for f in [18]:
    s=str(f)
    data=io.loadmat('E:/workshop/data/'+s+'/pcbf.mat')
    data=data['all_data_list']
    data=np.real(data)
    indices=io.loadmat('E:/workshop/data/'+s+'/pcbf_indices.mat')
    indices=indices['all_data_index']
    imus=io.loadmat('E:/workshop/data/'+s+'/imu_data.mat')
    imus=imus['imu_data']
    odmo=io.loadmat('E:/workshop/data/'+s+'/odom_data.mat')
    odmo=odmo['odom_data']

    x_last = 0
    p_last = 0
    Q = 0.5  #系统噪声
    Rk = 0.5  #测量噪声

    xs=[0];ys=[0];zs=[0];vxs=[0];vys=[0];vzs=[0];num_candidates=[0];vxx=[];vyy=[]
    x=0;y=0;z=0
    thr=0.1
    x_vec=np.array([1,0,0]);y_vec=np.array([0,1,0]);z_vec=np.array([0,0,1])
    num=140
    d_thr=0.15
    false_indices=[]
    for i in range(len(indices)-1):
        cur_data=data[indices[i,0]:indices[i+1,0],:]
        m=3;vx=0;vy=0;vz=0
        alphas=[];betas=[];gammas=[]
        candidates=[];real_candidates=[]
        if len(cur_data)>3:
            for k in range(len(cur_data)):
                    veck=np.array(cur_data[k,0:3])/math.sqrt(cur_data[k,0]**2+cur_data[k,1]**2+cur_data[k,2]**2)
                    theta_kx=math.acos(np.dot(veck,x_vec))
                    theta_ky=math.acos(np.dot(veck,y_vec))
                    theta_kz=math.acos(np.dot(veck,z_vec))
                    alphas.append(math.cos(theta_kx))
                    betas.append(math.cos(theta_ky))
                    gammas.append(math.cos(theta_kz))
                    
            for j in range(100):
                random.seed(j)
                aa,bb,cc = random.sample([o for o in range(len(cur_data))], 3)

                veca=np.array(cur_data[aa,0:3])/math.sqrt(cur_data[aa,0]**2+cur_data[aa,1]**2+cur_data[aa,2]**2)
                vecb=np.array(cur_data[bb,0:3])/math.sqrt(cur_data[bb,0]**2+cur_data[bb,1]**2+cur_data[bb,2]**2)
                vecc=np.array(cur_data[cc,0:3])/math.sqrt(cur_data[cc,0]**2+cur_data[cc,1]**2+cur_data[cc,2]**2)
                
                theta_ax=math.acos(np.dot(veca,x_vec))
                theta_ay=math.acos(np.dot(veca,y_vec))
                theta_az=math.acos(np.dot(veca,z_vec))
                theta_bx=math.acos(np.dot(vecb,x_vec))
                theta_by=math.acos(np.dot(vecb,y_vec))
                theta_bz=math.acos(np.dot(vecb,z_vec))
                theta_cx=math.acos(np.dot(vecc,x_vec))
                theta_cy=math.acos(np.dot(vecc,y_vec))
                theta_cz=math.acos(np.dot(vecc,z_vec))
                if theta_ax==theta_bx or theta_ax==theta_cx or theta_bx==theta_cx or theta_ay==theta_by or theta_ay==theta_cy or theta_by==theta_cy:
                    continue
                A=np.matrix([[math.cos(theta_ax),math.cos(theta_ay),math.cos(theta_az)],[math.cos(theta_bx),math.cos(theta_by),math.cos(theta_bz)],[math.cos(theta_cx),math.cos(theta_cy),math.cos(theta_cz)]])
                VR=np.matrix([[-cur_data[aa,3]],[-cur_data[bb,3]],[-cur_data[cc,3]]])
                [vxt,vyt,vzt]=A.I*VR
                vxt=float(vxt);vyt=float(vyt);vzt=float(vzt)
                dis_list=[cur_data[aa,0:3],cur_data[bb,0:3],cur_data[cc,0:3]]

                a2=alphas[aa]**2+alphas[bb]**2+alphas[cc]**2
                b2=betas[aa]**2+betas[bb]**2+betas[cc]**2
                c2=gammas[aa]**2+gammas[bb]**2+gammas[cc]**2
                ac=alphas[aa]*gammas[aa]+alphas[bb]*gammas[bb]+alphas[cc]*gammas[cc]
                bc=betas[aa]*gammas[aa]+betas[bb]*gammas[bb]+betas[cc]*gammas[cc]
                ab=alphas[aa]*betas[aa]+alphas[bb]*betas[bb]+alphas[cc]*betas[cc]
                ar=alphas[aa]*-cur_data[aa,3]+alphas[bb]*-cur_data[bb,3]+alphas[cc]*-cur_data[cc,3]
                br=betas[aa]*-cur_data[aa,3]+betas[bb]*-cur_data[bb,3]+betas[cc]*-cur_data[cc,3]
                cr=gammas[aa]*-cur_data[aa,3]+gammas[bb]*-cur_data[bb,3]+gammas[cc]*-cur_data[cc,3]
                
                dis_list=[cur_data[aa,0:3],cur_data[bb,0:3],cur_data[cc,0:3]]

                vote=0;candidates=[]
                for k in range(len(cur_data)):
                    if cur_data[k,1]>0.2:
                        v_cal=vxt*alphas[k]+vyt*betas[k]+vzt*gammas[k]
                        tmp_res=abs(v_cal+cur_data[k,3])
                        if tmp_res<thr:
                            test_dis=copy.deepcopy(dis_list)
                            test_dis=test_dis-cur_data[k,0:3]
                            flag=True
                            for o in range(len(test_dis[:,0])):
                                d=np.linalg.norm(test_dis[o,:])
                                if d<d_thr:
                                    flag=False
                            
                            if flag:
                                vote+=1
                                dis_list.append(cur_data[k,0:3])
                                dis_list.append(cur_data[k,0:3])
                                candidates.append(k)
                                a2=a2+alphas[k]**2
                                b2=b2+betas[k]**2
                                c2=c2+gammas[k]**2
                                ab=ab+alphas[k]*betas[k]
                                bc=bc+gammas[k]*betas[k]
                                ac=ac+gammas[k]*alphas[k]
                                ar=ar+alphas[k]*-cur_data[k,3]
                                br=br+betas[k]*-cur_data[k,3]
                                cr=cr+gammas[k]*-cur_data[k,3]

                if vote>m:  
                    real_candidates=copy.deepcopy(candidates)                
                    m=vote
                    L=np.matrix([[a2,ab,ac],[ab,b2,bc],[ac,bc,c2]])
                    R=np.matrix([[ar],[br],[cr]])
                    [vx,vy,vz]=L.I*R
                    vx=float(vx);vy=float(vy);vz=float(vz)
        if m<4:
            vxs.append(vxs[-1]); vys.append(vys[-1]); vzs.append(vzs[-1])
            false_indices.append(i)
        else:
            vxs.append(vx); vys.append(vy); vzs.append(vz)
        num_candidates.append(len(real_candidates))
        
        print("HERE",i,vx,vy,vz,len(cur_data),real_candidates)


    vz_tmp=[]
    for i in range(len(vzs)):
        pred,p_last,x_last = kalman(vzs[i],x_last,p_last,Q,Rk)
        vz_tmp.append(pred)
    vzs=vz_tmp  


    print("Length",len(vxs))


    dxs=[];dys=[];dzs=[]
    rotations=[]
    #xyz xzy yxz yzx zxy zyx
    for k in range(300):
        delta_x=sum(imus[0:(k+1)*20-1,3])*0.1/20
        delta_y=sum(imus[0:(k+1)*20-1,4])*0.1/20
        delta_z=sum(imus[0:(k+1)*20-1,5])*0.1/20
        rotations.append([delta_x,delta_y,delta_z])

    R0=np.matrix([[1,0,0],[0,1,0],[0,0,1]])

    vxt=[0];vyt=[0];vzt=[0]
    for k in range(300):
        ax=rotations[k][0];ay=rotations[k][1];az=rotations[k][2]
        Rx=np.matrix([[1,0,0],[0, math.cos(ax), -math.sin(ax)],[0, math.sin(ax), math.cos(ax)]])
        Ry=np.matrix([[math.cos(ay),0,math.sin(ay)],[0, 1, 0],[-math.sin(ay),0, math.cos(ay)]])
        Rz=np.matrix([[math.cos(az),-math.sin(az),0],[math.sin(az), math.cos(az),0],[0,0,1]])
        R0=-Rz*Ry*Rx
        vs= R0*np.matrix([[vys[k]],[vxs[k]],[vzs[k]]])
        vxt.append(float(vs[0]));vyt.append(float(vs[1]));vzt.append(float(vs[2]))
        x=x+float(vs[0])*0.1;y=y+float(vs[1])*0.1;z=z+float(vs[2])*0.1
        xs.append(x)
        ys.append(y)
        zs.append(z)
    xs=np.array(xs);ys=np.array(ys);zs=np.array(zs)

    # cs=[i for i in range(len(xs))]

    # data_index=f
    # pre = np.load(f'E:/drawing/output_array2/cross-11_9_1/{data_index}/pre_ep200.npy')
    # gt = np.load(f'E:/drawing/output_array2/cross-11_9_1/{data_index}/gt_ep200.npy')

    # errors_m=[]
    # fig=plt.figure()
    # for i in range(len(indices)-10):
    #     gts=gt[i+1, :] - gt[i, :]
    #     pres=pre[i+1, :] - pre[i, :]
    #     errors_m.append(np.linalg.norm(gts[0:2]-pres[0:2]))


    # errors=[];nums=[]
    # for i in range(len(indices)-2):
    #     gts=odmo[i+1,0:2] - odmo[i,0:2]
    #     pres=np.array([xs[i+1]-xs[i],ys[i+1]-ys[i]])
    #     errors.append(np.linalg.norm(gts-pres))
    #     nums.append(indices[i+1,0]-indices[i,0])
    # idx=[i for i in range(len(indices)-2)]
    # plt.subplot(2,1,1)
    # plt.title("Error")
    # plt.plot(idx,errors,label='Direct')
    # plt.plot([i for i in range(len(errors_m))],errors_m,label='milliego')
    # plt.legend()
    # plt.subplot(2,1,2)
    # # plt.title("Num of Points")
    # plt.plot(idx,num_candidates[0:299],label=str(data_index))
    # plt.legend()
    # plt.xlabel('Num of Points')
    # plt.show()
    # cm = plt.get_cmap('Reds')
    # cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    # scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    for i in range(100):
        fig = plt.figure(i)
        # plt.scatter(xs[false_indices],zs[false_indices],s=5,c='g')
        plt.plot(-odmo[0:(i*3)-1,0],odmo[0:(i*3)-1,1],linewidth= 4,label='Ground Truth')
        plt.plot(-xs[0:(i*3)-1],ys[0:(i*3)-1],linewidth= 4,label='Our Method')
        plt.legend(loc='upper right',prop=font2)
        plt.xlim([-6,0])
        plt.ylim([-1,5])
        # plt.xlabel('x')
        # plt.ylabel('y')
        plt.savefig(str(i)+'ransac_xy.png')
        # plt.show()

    # fig = plt.figure()
    # plt.scatter(xs[false_indices],zs[false_indices],s=5,c='g')
    # plt.plot(odmo[:,0],odmo[:,2])
    # plt.plot(xs,zs)
    # plt.xlabel('x')
    # plt.ylabel('z')
    # plt.savefig(s+'ransac_xz.jpg')

    # plt.show()

    # ax = Axes3D(fig)
    # ax.scatter(xs,ys,zs,s=5,c='b')
    # plt.scatter([0],[0],[0],c='r')
    # ax.plot(xs,ys,zs)
    # ax.set_xlim([-50,50])
    # ax.set_ylim([-50,50])
    # ax.set_zlim([-50,50])
    # # 添加坐标轴(顺序是Z, Y, X)


    # ax.set_zlabel('Z', fontdict={'size': 15, 'color': 'red'})
    # ax.set_ylabel('Y', fontdict={'size': 15, 'color': 'red'})
    # ax.set_xlabel('X', fontdict={'size': 15, 'color': 'red'})
    # plt.scatter(xs,zs,s=5)
    # ax.plot(xs[false_indices],ys[false_indices],zs[false_indices])
    

    # ovx=[0];ovy=[0];ovz=[0]
    # ovx.append((odmo[0][0]-0))
    # ovy.append((odmo[0][1]-0))
    # ovz.append((odmo[0][2]-0))
    # for i in range(1,300):
    #     ovx.append((odmo[i][0]-odmo[i-1][0]))
    #     ovy.append((odmo[i][1]-odmo[i-1][1]))
    #     ovz.append(-(odmo[i][2]-odmo[i-1][2]))
    # ovx=np.array(ovx);ovx=np.array(ovx);ovx=np.array(ovx);vxt=np.array(vxt)*0.1;vyt=np.array(vyt)*0.1;vzt=np.array(vzt)*0.1
    # plt.subplot(2,1,1)
    # plt.plot([i for i in range(301)],ovx)
    # plt.plot([i for i in range(301)],vxt)
    # plt.xlabel('Time')
    # plt.ylabel('Vx')
    # plt.subplot(2,1,2)
    # plt.plot([i for i in range(301)],abs(ovx-vxt))
    # plt.ylim([0,0.3])
    # plt.show()
    # plt.plot([i for i in range(301)],ovy)
    # plt.plot([i for i in range(301)],vyt)
    # plt.xlabel('Time')
    # plt.ylabel('Vy')
    # plt.show()
    # plt.plot([i for i in range(301)],ovz)
    # plt.plot([i for i in range(301)],vzt)
    # plt.xlabel('Time')
    # plt.ylabel('Vz')
    # plt.show()
    # for i in range(len(vxs)):
    #     if abs(vxs[i])>5 or abs(vys[i])>5 or abs(vzs[i])>5:
    #         print(vxs[i],vys[i],vzs[i],num_candidates[i])

    # print(np.mean(num_candidates))
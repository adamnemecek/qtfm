q = [
     quaternion(0.9584282037129221, 0.7300283992775096, 0.6470752153849828, 0.6470752153849828),
     quaternion(0.3877889209764178, 0.002728533663862054, 0.059233118245990646, 0.059233118245990646),
     quaternion(0.7095154552908696, 0.8519516087169816, 0.8613032947358639, 0.8613032947358639),
     quaternion(0.5527711909258979, 0.8388344201868314, 0.5105040217675825, 0.5105040217675825),
     quaternion(0.8561529554768845, 0.08150113997933417, 0.12941979046216212, 0.12941979046216212),
     quaternion(0.3006489646813413, 0.973174457685066, 0.9394428396213388, 0.9394428396213388),
     quaternion(0.886347907337274, 0.938496029682496, 0.4692614922066758, 0.4692614922066758),
     quaternion(0.7646076058356179, 0.4210837645964016, 0.6172734655044521, 0.6172734655044521),
     quaternion(0.7845772142254861, 0.9288234257080297, 0.5408316515442586, 0.5408316515442586),
     quaternion(0.5852633232268152, 0.4866818734898606, 0.6388974168100987, 0.6388974168100987)
 ]
 
 

basis = quaternion(0,1,1,1) / sqrt(3)



% a = [1,2,3]
default = unit(quaternion(1,1,1))

res = qdft2(q, default, 'L')

show(res)
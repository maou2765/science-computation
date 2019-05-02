def image_dist(f,do):
    return f*do/(do-f)
def object_dist_v(f,do,separation):
    return separation-image_dist(f,do)
def magnification(f,do):
    return f/(f-do)
dist=[object_dist_v(-22,50,104)]
dist.append(object_dist_v(48.,dist[-1],110))
dist.append(object_dist_v(48.,dist[-1],686-491))
print(-1*magnification(-22,50.))
print(-1*magnification(-22,50.)*magnification(48.,dist[0]))
print(-1*magnification(-22,50.)*magnification(48.,dist[0])*magnification(48,dist[1]))
print(-1*magnification(-22,50.)*magnification(48.,dist[0])*magnification(48,dist[1])*magnification(48,dist[2]))

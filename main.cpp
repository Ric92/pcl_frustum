#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>
#include <pcl/filters/frustum_culling.h>

class Sample
{
public:
  static int uniform(int from, int to);
  static double uniform();
  static double gaussian(double sigma);
};

static double uniform_rand(double lowerBndr, double upperBndr)
{
  return lowerBndr + ((double)std::rand() / (RAND_MAX + 1.0)) * (upperBndr - lowerBndr);
}

static double gauss_rand(double mean, double sigma)
{
  double x, y, r2;
  do
  {
    x = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
    y = -1.0 + 2.0 * uniform_rand(0.0, 1.0);
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0.0);
  return mean + sigma * y * std::sqrt(-2.0 * log(r2) / r2);
}

int Sample::uniform(int from, int to)
{
  return static_cast<int>(uniform_rand(from, to));
}

double Sample::uniform()
{
  return uniform_rand(0., 1.);
}

double Sample::gaussian(double sigma)
{
  return gauss_rand(0., sigma);
}

void drawFrustum(pcl::FrustumCulling<pcl::PointXYZRGB> _fc, boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, int id)
{

  Eigen::Matrix4f cameraPose = _fc.getCameraPose();
  float hfov_ = _fc.getHorizontalFOV();
  float vfov_ = _fc.getVerticalFOV();
  float np_dist_ = _fc.getNearPlaneDistance();
  float fp_dist_ = _fc.getFarPlaneDistance();

  Eigen::Vector4f pl_n; // near plane
  Eigen::Vector4f pl_f; // far plane
  Eigen::Vector4f pl_t; // top plane
  Eigen::Vector4f pl_b; // bottom plane
  Eigen::Vector4f pl_r; // right plane
  Eigen::Vector4f pl_l; // left plane

  Eigen::Vector3f view = cameraPose.block(0, 0, 3, 1);  // view vector for the camera  - first column of the rotation matrix
  Eigen::Vector3f up = cameraPose.block(0, 1, 3, 1);    // up vector for the camera    - second column of the rotation matrix
  Eigen::Vector3f right = cameraPose.block(0, 2, 3, 1); // right vector for the camera - third column of the rotation matrix
  Eigen::Vector3f T = cameraPose.block(0, 3, 3, 1);     // The (X, Y, Z) position of the camera w.r.t origin
  float vfov_rad = float(vfov_ * M_PI / 180);           // degrees to radians
  float hfov_rad = float(hfov_ * M_PI / 180);           // degrees to radians

  float np_h = float(2 * tan(vfov_rad / 2) * np_dist_); // near plane height
  float np_w = float(2 * tan(hfov_rad / 2) * np_dist_); // near plane width

  float fp_h = float(2 * tan(vfov_rad / 2) * fp_dist_); // far plane height
  float fp_w = float(2 * tan(hfov_rad / 2) * fp_dist_); // far plane width

  Eigen::Vector3f fp_c(T + view * fp_dist_);                          // far plane center
  Eigen::Vector3f fp_tl(fp_c + (up * fp_h / 2) - (right * fp_w / 2)); // Top left corner of the far plane
  Eigen::Vector3f fp_tr(fp_c + (up * fp_h / 2) + (right * fp_w / 2)); // Top right corner of the far plane
  Eigen::Vector3f fp_bl(fp_c - (up * fp_h / 2) - (right * fp_w / 2)); // Bottom left corner of the far plane
  Eigen::Vector3f fp_br(fp_c - (up * fp_h / 2) + (right * fp_w / 2)); // Bottom right corner of the far plane

  Eigen::Vector3f np_c(T + view * np_dist_);                          // near plane center
  Eigen::Vector3f np_tl(np_c + (up * np_h / 2) - (right * np_w / 2)); // Top left corner of the near plane
  Eigen::Vector3f np_tr(np_c + (up * np_h / 2) + (right * np_w / 2)); // Top right corner of the near plane
  Eigen::Vector3f np_bl(np_c - (up * np_h / 2) - (right * np_w / 2)); // Bottom left corner of the near plane
  Eigen::Vector3f np_br(np_c - (up * np_h / 2) + (right * np_w / 2)); // Bottom right corner of the near plane

  pl_f.block(0, 0, 3, 1).matrix() = (fp_bl - fp_br).cross(fp_tr - fp_br); // Far plane equation - cross product of the
  pl_f(3) = -fp_c.dot(pl_f.block(0, 0, 3, 1));                            // perpendicular edges of the far plane

  pl_n.block(0, 0, 3, 1).matrix() = (np_tr - np_br).cross(np_bl - np_br); // Near plane equation - cross product of the
  pl_n(3) = -np_c.dot(pl_n.block(0, 0, 3, 1));                            // perpendicular edges of the far plane

  Eigen::Vector3f a(fp_bl - T); // Vector connecting the camera and far plane bottom left
  Eigen::Vector3f b(fp_br - T); // Vector connecting the camera and far plane bottom right
  Eigen::Vector3f c(fp_tr - T); // Vector connecting the camera and far plane top right
  Eigen::Vector3f d(fp_tl - T); // Vector connecting the camera and far plane top left

  pcl::PointXYZRGB camera;
  camera.x = T[0];camera.y = T[1];camera.z = T[2];

  pcl::PointXYZRGB ntl; ntl.x = np_tl[0]; ntl.y = np_tl[1]; ntl.z = np_tl[2];
  pcl::PointXYZRGB ntr; ntr.x = np_tr[0]; ntr.y = np_tr[1]; ntr.z = np_tr[2];
  pcl::PointXYZRGB nbl; nbl.x = np_bl[0]; nbl.y = np_bl[1]; nbl.z = np_bl[2];
  pcl::PointXYZRGB nbr; nbr.x = np_br[0]; nbr.y = np_br[1]; nbr.z = np_br[2];
  pcl::PointXYZRGB ftl; ftl.x = fp_tl[0]; ftl.y = fp_tl[1]; ftl.z = fp_tl[2];
  pcl::PointXYZRGB ftr; ftr.x = fp_tr[0]; ftr.y = fp_tr[1]; ftr.z = fp_tr[2];
  pcl::PointXYZRGB fbl; fbl.x = fp_bl[0]; fbl.y = fp_bl[1]; fbl.z = fp_bl[2];
  pcl::PointXYZRGB fbr; fbr.x = fp_br[0]; fbr.y = fp_br[1]; fbr.z = fp_br[2];

  viewer->addCoordinateSystem(1.0, T[0], T[1], T[2], "camera" + std::to_string(id));

  viewer->addLine(camera, ftl, 0.0, 0.0, 1.0, "camera_fartopleft" + std::to_string(id));
  viewer->addLine(camera, ftr, 0.0, 0.0, 1.0, "camera_fartopright" + std::to_string(id));
  viewer->addLine(camera, fbl, 1.0, 0.0, 0.0, "camera_nearbotleft" + std::to_string(id));
  viewer->addLine(camera, fbr, 1.0, 0.0, 0.0, "camera_nearbotright" + std::to_string(id));

  viewer->addLine(ftl, ftr, 1.0, 0.0, 0.0, "fartopleft_fartopright" + std::to_string(id));
  viewer->addLine(ftl, fbl, 1.0, 0.0, 0.0, "fartopleft_farbotleft" + std::to_string(id));
  viewer->addLine(fbr, ftr, 1.0, 0.0, 0.0, "farbotright_fartopright" + std::to_string(id));
  viewer->addLine(fbr, fbl, 1.0, 0.0, 0.0, "farbotright_farbotleft" + std::to_string(id));

  viewer->addLine(ntl, ntr, 1.0, 0.0, 0.0, "neartopleft_neartopright" + std::to_string(id));
  viewer->addLine(ntl, nbl, 1.0, 0.0, 0.0, "neartopleft_nearbopleft" + std::to_string(id));
  viewer->addLine(nbr, ntr, 1.0, 0.0, 0.0, "nearbotright_neartopright" + std::to_string(id));
  viewer->addLine(nbr, nbl, 1.0, 0.0, 0.0, "nearbotright_nearbopleft" + std::to_string(id));
}

int main(int argc, char **argv)
{
  boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
  viewer->setBackgroundColor(100, 100, 100);
  viewer->addCoordinateSystem(1.0, 0.0, 0.0, 0.0, "camera");
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);

  std::vector<Eigen::Vector3d> true_points;
  for (size_t i = 0; i < 500; ++i)
  {
    Eigen::Vector3d point;
    point = {(float)((Sample::uniform() + 0.5) * 3),
             ((float)(Sample::uniform() - 0.5) * 3),
             ((float)(Sample::uniform() - 0.5) * 3)};

    true_points.push_back(Eigen::Vector3d((double)point[0],
                                          (double)point[1],
                                          (double)point[2]));

    pcl::PointXYZRGB p(0, 0, 255);
    p.x = point[0];
    p.y = point[1];
    p.z = point[2];

    cloud->push_back(p);
  }

  viewer->addPointCloud<pcl::PointXYZRGB>(cloud, "firstCloud");
  viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "firstCloud");
  viewer->spinOnce(10, true);

  pcl::FrustumCulling<pcl::PointXYZRGB> fc;
  fc.setInputCloud(cloud);
  fc.setVerticalFOV(60);
  fc.setHorizontalFOV(90);
  fc.setNearPlaneDistance(1.0);
  fc.setFarPlaneDistance(3.0);
  Eigen::Matrix4f camera_pose;
  camera_pose.setIdentity();
  fc.setCameraPose(camera_pose);
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr target(new pcl::PointCloud<pcl::PointXYZRGB>);
  fc.filter(*target);
  viewer->addPointCloud<pcl::PointXYZRGB>(target, "secondCloud");
  viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "secondCloud");
  viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "secondCloud");

  drawFrustum(fc, viewer,1);

  std::cout << "Filtered " + std::to_string(target->size()) + " points of " + std::to_string(cloud->size()) << std::endl;

  pcl::FrustumCulling<pcl::PointXYZRGB> fc2;
  fc.setInputCloud(target);
  fc.setVerticalFOV(60);
  fc.setHorizontalFOV(90);
  fc.setNearPlaneDistance(1.0);
  fc.setFarPlaneDistance(3.0);
  Eigen::Matrix4f camera2_pose;
  camera2_pose.setIdentity();
  camera2_pose.col(3).head<3>() << 1,1,1;
  fc.setCameraPose(camera2_pose);
  pcl::PointCloud<pcl::PointXYZRGB>::Ptr target2(new pcl::PointCloud<pcl::PointXYZRGB>);
  fc.filter(*target2);
  viewer->addPointCloud<pcl::PointXYZRGB>(target2, "thirdCloud");
  viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "thirdCloud");
  viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "thirdCloud");

  drawFrustum(fc, viewer,2);

  std::cout << "Filtered " + std::to_string(target2->size()) + " points of " + std::to_string(target->size()) << std::endl;

  int a = -3;
  while (a < 0)
  {
    viewer->spinOnce(10, true);
  }

  return (0);
}

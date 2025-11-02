#include <iostream>
#include <gdal_priv.h>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_tif_file>" << std::endl;
        return 1;
    }

    const char* filepath = argv[1];
    std::cout << "--- GDAL Minimal Test ---" << std::endl;

    // 1. 测试GDAL初始化
    std::cout << "Step 1: Calling GDALAllRegister()..." << std::endl;
    GDALAllRegister();
    std::cout << "Step 1: OK. GDAL initialized." << std::endl;

    // 2. 测试文件打开
    std::cout << "Step 2: Calling GDALOpen() with path: " << filepath << "..." << std::endl;
    GDALDataset* dataset = (GDALDataset*)GDALOpen(filepath, GA_ReadOnly);

    if (dataset == nullptr) {
        std::cout << "Step 2: FAILED. GDALOpen returned nullptr." << std::endl;
        return 1;
    }
    std::cout << "Step 2: OK. File opened successfully." << std::endl;

    // 3. 测试读取元数据
    std::cout << "Step 3: Reading raster size..." << std::endl;
    int nXSize = dataset->GetRasterXSize();
    int nYSize = dataset->GetRasterYSize();
    std::cout << "Step 3: OK. Size is " << nXSize << " x " << nYSize << "." << std::endl;

    // 4. 清理
    std::cout << "Step 4: Closing dataset..." << std::endl;
    GDALClose(dataset);
    std::cout << "Step 4: OK. Test finished successfully." << std::endl;

    return 0;
}
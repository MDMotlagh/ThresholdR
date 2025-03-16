The bmcite dataset is a CITE-seq experiment on bone marrow samples. It is available through the SeuratData R package. 
```r
SeuratData::AvailableData()
SeuratData::InstallData('bmcite')

#then load the data
so <- SeuratData::LoadData('bmcite')
```

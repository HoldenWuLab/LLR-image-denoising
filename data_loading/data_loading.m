function [kdata, acs, prescan] = data_loading(data_type)

if data_type == 1
    load("./example_dataset/phantom/acs.mat","acs");
    load("./example_dataset/phantom/prescan.mat","prescan");
    load("./example_dataset/phantom/kdata_ec1.mat","kdata_ec1");
    load("./example_dataset/phantom/kdata_ec2.mat","kdata_ec2");
    load("./example_dataset/phantom/kdata_ec3.mat","kdata_ec3");
    load("./example_dataset/phantom/kdata_ec4.mat","kdata_ec4");
    load("./example_dataset/phantom/kdata_ec5.mat","kdata_ec5");
    load("./example_dataset/phantom/kdata_ec6.mat","kdata_ec6");
    kdata = cat(5,kdata_ec1,kdata_ec2,kdata_ec3,kdata_ec4,kdata_ec5,kdata_ec6);
elseif data_type == 2
    load("./example_dataset/liver/acs.mat","acs");
    load("./example_dataset/liver/prescan.mat","prescan");
    load("./example_dataset/liver/kdata_ec1.mat","kdata_ec1");
    load("./example_dataset/liver/kdata_ec2.mat","kdata_ec2");
    load("./example_dataset/liver/kdata_ec3.mat","kdata_ec3");
    load("./example_dataset/liver/kdata_ec4.mat","kdata_ec4");
    load("./example_dataset/liver/kdata_ec5.mat","kdata_ec5");
    load("./example_dataset/liver/kdata_ec6.mat","kdata_ec6");
    kdata = cat(5,kdata_ec1,kdata_ec2,kdata_ec3,kdata_ec4,kdata_ec5,kdata_ec6);
end
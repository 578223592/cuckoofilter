#include "cuckoofilter.h"

#include <assert.h>
#include <math.h>

#include <chrono>
#include <iostream>
#include <vector>

using cuckoofilter::CuckooFilter;




int main(int argc, char **argv) {
  size_t total_items = 1000000;
  // size_t total_items = 100;

  // Create a cuckoo filter where each item is of type size_t and
  // use 12 bits for each item:
  //    CuckooFilter<size_t, 12> filter(total_items);
  // To enable semi-sorting, define the storage of cuckoo filter to be
  // PackedTable, accepting keys of size_t type and making 13 bits
  // for each key:
  //   CuckooFilter<size_t, 13, cuckoofilter::PackedTable> filter(total_items);
  CuckooFilter<size_t, 12> filter(total_items);


  // Insert items to this cuckoo filter
  size_t num_inserted = 0;
  for (size_t i = 0; i < total_items; i++, num_inserted++) {
    if (filter.Add(i) != cuckoofilter::Ok) {
      std::cerr << "Insertion failed at " << i << std::endl;
      break;
    }
  }

  // Check if previously inserted items are in the filter, expected
  // true for all items
  for (size_t i = 0; i < num_inserted; i++) {
    assert(filter.Contain(i) == cuckoofilter::Ok);
  }

  // Check non-existing items, a few false positives expected
  // 通过插入全部不存在item来检测FP率
  size_t total_queries = 0;
  size_t false_queries = 0;
  for (size_t i = total_items; i < 2 * total_items; i++) {
    if (filter.Contain(i) == cuckoofilter::Ok) {
      false_queries++;
    }
    total_queries++;
  }


  // Output the measured false positive rate
  std::cout << "false positive rate is "
            << 100.0 * false_queries / total_queries << "%\n";

  std::cout<<filter.Info()<<std::endl;

  std::cout<<"占用空间大小："<<sizeof(filter)<<std::endl;



  // hash 和 截取长度的速度测试
  CuckooFilter<std::string, 12> filter2(total_items);
  size_t yyy = 0;
  std::string y = "12345678901234567890";
  size_t index;uint32_t tag;


  int i = 0;
  // 开始计时
  auto start = std::chrono::high_resolution_clock::now();


  for(;i<5000;++i) {
    // cuckoofilter::HashUtil::MurmurHash(y);
    cuckoofilter::HashUtil::SuperFastHash(y);
    // std::string tmp = y.substr(0,12);
  }
  // filter.GenerateIndexTagHashTest(yyy, &index, &tag);
  // cuckoofilter::HashUtil::MurmurHash(y);
  // filter2.GenerateIndexTagHashTest(y,&index,&tag);

  // 结束计时
  auto end = std::chrono::high_resolution_clock::now();

  // 计算耗时并输出结果
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
  std::cout << "Function execution time: " << duration.count()/ (i+1)<< " nanoseconds." << std::endl;

  return 0;
}

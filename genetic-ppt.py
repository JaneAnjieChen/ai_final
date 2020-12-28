# original version
import random
import numpy as np
import copy
import matplotlib.pyplot as plt

# 初始化种群，赋予每个个体一个全局唯一的编号
def ini_species(gene_num, ini_size):
    population = {}
    id = 0
    for i in range(ini_size):
        id += 1
        individual = []
        individual_info = {}
        for gene in range(gene_num):
            individual.append(random.randint(0, 1))

        individual_info['gene seq'] = individual
        population[id] = individual_info
    # print('initial population: ', population)
    return population


# 设定个体的适应度函数f(x)
def fitness_f(x):
    return (15*x - x*x)


# 将一个0、1列表转化成一个二进制字符串
def get_HEX(l):
    hex_str = ''
    for i in l:
        hex_str += str(i)
    return hex_str


# 计算一个种群字典中所有个体的适应度值，计算每个个体被选择的概率Ps,
# 并返回扩充了适应度值、选择概率后的种群字典
def calc_fitness_Ps(population):
    total_fitness = 0
    for id, individual_info in population.items():
        gene_dec = int(get_HEX(individual_info['gene seq']), 2)
        individual_info['fitness'] = fitness_f(gene_dec)
        total_fitness += fitness_f(gene_dec)
    # print(population)

    for id, individual_info in population.items():
        individual_info['Ps'] = individual_info['fitness']/total_fitness
    # print(population)

    return population


# 传入参数：{1: {'gene seq': , 'fitness': , 'Ps': }}
# Selection，轮盘赌随机算法
def selection(population):
    cumulation_Ps = 0
    for id, individual_info in population.items():
        individual_info['cumulated Ps'] = cumulation_Ps
        cumulation_Ps += individual_info['Ps']
        individual_info['selected times'] = 0
    # print(population)

    for i in range(len(population)):
        rand_float = random.random()
        for id, individual_info in population.items():
            if rand_float < individual_info['cumulated Ps']:
                # found ceiling，这个节点的前一个节点便是下界
                population[id-1]['selected times'] += 1
                break
    # print(population)

    # 找到了selection次数，现在要进行配对
    return population


# 传入参数：{id: {'selected times': }}
# return 配对列表
def matching(population):
    # [[], []]
    match_seq = []
    tmp_parents = []
    count = 1

    while count != 0:
        count = 0
        for id, individual_info in population.items():
            if len(tmp_parents) < 2 and individual_info['selected times'] != 0:
                tmp_parents.append(id)
                individual_info['selected times'] -= 1
                count += 1

            if len(tmp_parents) == 2:
                match_seq.append(tmp_parents)
                tmp_parents = []

    # print(match_seq)
    return match_seq


# 一对染色体的交叉：随机选择断裂点，交叉，返回得到的一对新染色体
# 传入参数：两个list
def crossover(gene_seq1, gene_seq2):
    if len(gene_seq1) != len(gene_seq2):
        print('这两个染色体不一样长，出错，程序退出')
        exit(0)

    # 断裂点索引
    break_index = random.randint(0, len(gene_seq1))

    gene_seq1_0 = gene_seq1[:break_index]
    gene_seq1_1 = gene_seq1[break_index:]
    gene_seq2_0 = gene_seq2[:break_index]
    gene_seq2_1 = gene_seq2[break_index:]

    new_gene_seq1 = gene_seq1_0 + gene_seq2_1
    new_gene_seq2 = gene_seq2_0 + gene_seq1_1

    return new_gene_seq1, new_gene_seq2


# 随机选择一个基因（1位）进行突变
# 传入参数：一个list
def mutation(gene_seq):
    mutation_index = random.randint(0, len(gene_seq)-1)
    if gene_seq[mutation_index] == 0:
        gene_seq[mutation_index] = 1
    else:
        gene_seq[mutation_index] = 0
    return gene_seq


# 计算平均适应度
def calc_ave_fitness(population):
    ave_fitness = 0
    total = 0
    count = 0
    for id, individual_info in population.items():
        total += individual_info['fitness']
        count += 1
    ave_fitness = total/count
    return ave_fitness

# def whether_found_optimal(population):
#     genes = []
#     for individual_id, individual_info in population.items():
#         genes.append(individual_info['gene seq'])


def propagation(gene_num, ini_size, crossover_p, mutation_p):
    # 重复繁殖过程，直到下一代种群数量达到N， 用新种群（其中均为新一代）替换初始种群。
    N = 300
    population = ini_species(gene_num, ini_size)
    found_optimal = False
    ave_fitness = []
    generation_count = 0
    while found_optimal == False:
        while len(population) <= N:
            population = calc_fitness_Ps(population)
            population = selection(population)
            ave_fitness.append(calc_ave_fitness(population))
            match_seq = matching(population)

            # 本次繁殖新生成的种群
            tmp_population = {}
            for parents in match_seq:
                # 按照crossover_p进行基因交叉
                if round(np.random.uniform(0, 1), 1) <= crossover_p:
                    new_gene_seq1, new_gene_seq2 = crossover(population[parents[0]]['gene seq'], population[parents[1]]['gene seq'])
                else:
                    new_gene_seq1, new_gene_seq2 = population[parents[0]]['gene seq'], population[parents[1]]['gene seq']

                # 按照mutation_p进行基因突变
                if round(np.random.uniform(0, 1), 1) <= mutation_p:
                    new_gene_seq1 = mutation(new_gene_seq1)
                if round(np.random.uniform(0, 1), 1) <= mutation_p:
                    new_gene_seq2 = mutation(new_gene_seq2)


                # 将新生成的个体加入到population字典中和tmp_population中以防种群大小达到N需要替换初始种群
                end_id = list(population.keys())[-1]
                population[end_id+1] = {}
                population[end_id+2] = {}
                population[end_id+1]['gene seq'] = new_gene_seq1
                population[end_id+2]['gene seq'] = new_gene_seq2

                tmp_population[end_id+1] = {}
                tmp_population[end_id+2] = {}
                tmp_population[end_id+1]['gene seq'] = new_gene_seq1
                tmp_population[end_id+2]['gene seq'] = new_gene_seq2
            generation_count += 1

        population = copy.deepcopy(tmp_population)

        # 收敛，有接近最优解，found_optimal置为True
        if generation_count >= 1000:
            found_optimal = True

    # 画图
    plt.plot(range(generation_count), ave_fitness, 'b')
    plt.xlabel('Generation')
    plt.ylabel('Average Fitness')
    plt.show()


    # population = ini_species(gene_num, ini_size)
    # generation_count = 0
    # find_optimal_flag = False
    # ave_fitness = []
    # while find_optimal_flag == False:
    #     population = calc_fitness_Ps(population)
    #     print(population)
    #     print('==================================================================')
    #
    #     generation_count += 1
    #     ave_fitness.append(calc_ave_fitness(population))
    #     # 重复繁殖过程，直到下一代种群数量达到N， 用新种群（其中均为新一代）替换初始种群。
    #     # 改成
    #     if generation_count == 20:
    #         find_optimal_flag = True
    #         break
    #
    #     # 选择次数定下来了
    #     population = selection(population)
    #     # 配对
    #     match_seq = matching(population)
    #
    #     for parents in match_seq:
    #         # 按照crossover_p进行基因交叉
    #         if round(np.random.uniform(0, 1), 1) <= crossover_p:
    #             new_gene_seq1, new_gene_seq2 = crossover(population[parents[0]]['gene seq'], population[parents[1]]['gene seq'])
    #         else:
    #             new_gene_seq1, new_gene_seq2 = population[parents[0]]['gene seq'], population[parents[1]]['gene seq']
    #
    #         # 按照mutation_p进行基因突变
    #         if round(np.random.uniform(0, 1), 1) <= mutation_p:
    #             new_gene_seq1 = mutation(new_gene_seq1)
    #         if round(np.random.uniform(0, 1), 1) <= mutation_p:
    #             new_gene_seq2 = mutation(new_gene_seq2)
    #
    #         # 将新生成的个体加入到population字典中
    #         end_id = list(population.keys())[-1]
    #         population[end_id+1] = {}
    #         population[end_id+2] = {}
    #         population[end_id+1]['gene seq'] = new_gene_seq1
    #         population[end_id+2]['gene seq'] = new_gene_seq2


if __name__ == '__main__':
    # 物种个体的基因总数
    # gene_num = int(input('The number of genes: '))
    gene_num = 4
    # 初始种群数量
    # ini_size = int(input('The size of the initial group: '))
    ini_size = 6

    # population = ini_species(gene_num, ini_size)
    #
    # population = calc_fitness_Ps(population)
    #
    # population = selection(population)
    #
    # matching(population)

    # 设定交叉概率Pc
    # crossover_p = float(input('Crossover probability: '))
    crossover_p = 0.7
    # 突变概率为Pm
    # mutation_p = float(input('Mutation rate: '))
    mutation_p = 0.001

    propagation(gene_num, ini_size, crossover_p, mutation_p)

### SIMULATION MAIN ALGORITHM

import numpy as np
from scipy.stats import erlang
from collections import deque
import matplotlib.pyplot as plt

# There are a few things left that are no longer needed
# parameters R_A, R_B, U_CPU, U_DISK (probably) don't provide any infomation and can be removed
# function status and whatever derives from its use is also leftover and doesn't provide anything useful


sec_of_msec = lambda t: t * 0.001

D_CPU_A: float = sec_of_msec(52.2911)
D_DISK_A: float = sec_of_msec(109.315)
v_DISK_A: float = 81.6255
D_OUT_A: float = sec_of_msec(8.84955)

D_CPU_B: float = sec_of_msec(4.41479)
D_DISK_B: float = sec_of_msec(140.627)
v_DISK_B: float = 56.978
D_OUT_B: float = sec_of_msec(78.6918)

D_CPU_C: float = sec_of_msec(59.6817)
D_DISK_C: float = sec_of_msec(164.246)
v_DISK_C: float = 69.7886
D_OUT_C: float = sec_of_msec(95.8399)

v_CPU_A: float = v_DISK_A + 1
v_CPU_B: float = v_DISK_B + 1
v_CPU_C: float = v_DISK_C + 1
v_OUT_A = v_OUT_B = v_OUT_C = 1

arrival_rate: float = 256776 / 45100  # λ: total jobs per total time frame
average_arrival_time: float = 1 / arrival_rate

station_index = lambda st: 0 if st == "CPU" else 1 if st == "DISK" else 2
category_index = lambda c: 0 if c == "A" else 1 if c == "B" else 2
Dij = [
    [D_CPU_A, D_CPU_B, D_CPU_C],
    [D_DISK_A, D_DISK_B, D_DISK_C],
    [D_OUT_A, D_OUT_B, D_OUT_C],
]
vij = [
    [v_CPU_A, v_CPU_B, v_CPU_C],
    [v_DISK_A, v_DISK_B, v_DISK_C],
    [v_OUT_A, v_OUT_B, v_OUT_C],
]
Sij = [[], [], []]
for i in range(3):
    for j in range(3):
        Sij[i].append(Dij[i][j] / vij[i][j])


def Erlang_4_service_time(category: str) -> float:
    # mean = k / λ
    k = 4
    lamda = 4 / Sij[station_index("CPU")][category_index(category)]
    return erlang.rvs(a=k, scale=1 / lamda)


def Exponential_service_time(station: str, category: str) -> float:
    # mean = 1 / λ
    mean = Sij[station_index(station)][category_index(category)]
    return np.random.exponential(scale=mean)


U = lambda: float(np.random.uniform())


def job_category() -> str:
    p_A = 0.30640
    p_B = 0.37634
    # p_C = 0.31726
    category_prob = U()
    if category_prob <= p_A:
        return "A"
    elif category_prob <= p_A + p_B:
        return "B"
    else:
        return "C"


# Initialize

# plot_list has 4 lists, one for each station. Every time an event happens, the number of jobs
# using the station, either waiting in its queue or not, is appended on the station's list.
# The last list is for the jobs that 'backed' from the system.
plot_list: list[list[int]] = [[], [], [], []]


def status(event: str) -> None:  # CHECK
    plot_list[0].append(len(STATION_queue[0]))
    for i in range(1, 3):  # i = i_DISK, i_OUT
        plot_list[i].append(0 if empty_STATION[i] else len(STATION_queue[i]) + 1)
    plot_list[3].append(backed_jobs)
    return


STATION_queue = [deque(), deque(), deque()]

STATION_current_job: list[int] = [None, None, None]
STATION_current_category: list[str] = [None, None, None]
STATION_start_time: list[float] = [None, None, None]
STATION_remaining_time: list[float] = [None, None, None]

STATION_EMPTY_TIME: list[float] = [0, 0, 0]
STATION_LAST_EMPTY: list[float] = [0, 0, 0]
empty_STATION: list[bool] = [True, True, True]


def set_current_counters(
    station: str, job: int, category: str, start_time: float, remaining_time: float
) -> None:
    i = station_index(station)
    STATION_current_job[i] = job
    STATION_current_category[i] = category
    STATION_start_time[i] = start_time
    STATION_remaining_time[i] = remaining_time
    return


valid_sum = lambda x, y: x + y if x is not None and y is not None else None
valid_sub = lambda x, y: x - y if x is not None and y is not None else None


def next_event(arrival_time: float) -> str:
    CPU_finish_time = valid_sum(STATION_start_time[0], STATION_remaining_time[0])
    DISK_finish_time = valid_sum(STATION_start_time[1], STATION_remaining_time[1])
    OUT_finish_time = valid_sum(STATION_start_time[2], STATION_remaining_time[2])

    lst = [arrival_time, CPU_finish_time, DISK_finish_time, OUT_finish_time]

    min_value = float("inf")
    min_position = None

    for index, value in enumerate(lst):
        if value is not None and value < min_value:
            min_value = value
            min_position = index

    return ["arrival", "CPU", "DISK", "OUT"][min_position]


# # arrival_time(job_id: int) -> arrival_time_of_job: float
arrival_time_of_job = lambda job_id: job_id * average_arrival_time
# departure_time: dict[job_id, departure_time_of_job]
departure_time: dict[int, float] = {}
# disk_visits: dict[job_id, current_number_of_visits]
disk_visits: dict[int, float] = {}

clock: int = 0
job_id: int = 1
curr_jobs: int = 0
backed_jobs: int = 0
theta = lambda: np.random.normal(loc=12, scale=3)

arrivals_dict = {}
arrivals_dict[average_arrival_time] = []
for i in range(3):
    for j in range(3):
        arrivals_dict[Dij[i][j]] = []


R_A_different = []
R_B_different = []
R_C_different = []
U_CPU_different = []
U_DISK_different = []
U_OUT_different = []


def decrease_CPU_remaining_time(time_passed: float) -> None:
    if not time_passed:
        return
    for job in STATION_queue[i]:
        job[2] -= time_passed
    return


def add_job_to_station(station: str, job_id: int, job_category: str) -> None:
    global clock
    i = station_index(station)
    service_time = float(
        Erlang_4_service_time(job_category)
        if station == "CPU"
        else Exponential_service_time(station, job_category)
    )

    def add_job_to_CPU_queue(
        job_id: int, job_category: str, service_time: float
    ) -> None:
        # nonlocal i
        i = 0  # CPU
        insert_position = len(STATION_queue[i])
        for index, value in enumerate(STATION_queue[i]):
            if value[2] > service_time:
                insert_position = index
                break
        STATION_queue[i].insert(insert_position, [job_id, job_category, service_time])
        return

    match station:
        case "CPU":
            decrease_CPU_remaining_time(valid_sub(clock, STATION_start_time[i]))
            add_job_to_CPU_queue(job_id, job_category, service_time)
            job, category, remaining_time = STATION_queue[i][0]
            set_current_counters(station, job, category, clock, remaining_time)
            if empty_STATION[i]:
                STATION_EMPTY_TIME[i] += clock - STATION_LAST_EMPTY[i]
                STATION_LAST_EMPTY[i] = clock
                empty_STATION[i] = False
        case _:
            if empty_STATION[i]:
                set_current_counters(station, job_id, job_category, clock, service_time)
                if empty_STATION[i]:
                    STATION_EMPTY_TIME[i] += clock - STATION_LAST_EMPTY[i]
                    STATION_LAST_EMPTY[i] = clock
                    empty_STATION[i] = False
            else:
                STATION_queue[i].append([job_id, job_category, service_time])

    return


def load_job_from_queue(station: str) -> None:
    global clock
    i = station_index(station)
    if not STATION_queue[i]:
        set_current_counters(station, None, None, None, None)
        empty_STATION[i] = True
        STATION_LAST_EMPTY[i] = clock
    else:
        job, category, remaining_time = STATION_queue[i].popleft()
        set_current_counters(station, job, category, clock, remaining_time)

    return


old_clock = 0
job_A, job_B, job_C = [], [], []

regens = 0
regen_cycle_length_list = []
num_of_A_jobs_in_regen_cycle = []
num_of_B_jobs_in_regen_cycle = []
num_of_C_jobs_in_regen_cycle = []
Checkpoint = True  # <10% μέσης τιμής

# average() computes the average value of a list disregarding its 1st element
average = lambda lst: sum(lst[1:]) / (len(lst) - 1)

while Checkpoint and regens < 1000:

    is_system_in_initial_state = empty_STATION == [True, True, True]
    if is_system_in_initial_state:
        # check the confidence interval every 20 cycles
        if (regens - 1) % 20 == 0 and regens != 1:

            # calculate s ^ 2
            y_bar = average(R_B_different)
            c_bar = average(num_of_B_jobs_in_regen_cycle)
            n = regens - 1

            sy = 0
            sc = 0
            syc = 0
            for i in range(1, regens):
                sy += (R_B_different[i] - y_bar) ** 2
                sc += (num_of_B_jobs_in_regen_cycle[i] - c_bar) ** 2
                syc += (R_B_different[i] - y_bar) * (
                    num_of_B_jobs_in_regen_cycle[i] - c_bar
                )

            sy /= n - 1
            sc /= n - 1
            syc /= n - 1
            R = average(R_B_different) / average(num_of_B_jobs_in_regen_cycle)
            s = sy - 2 * R * syc + R**2 * sc

            # z for 0.95
            z1_a_2 = 1.960
            confidence_interval = (
                z1_a_2 * float(np.sqrt(s)) / (c_bar * float(np.sqrt(n)))
            )

            if confidence_interval / R < 0.1:  # 10%
                Checkpoint = False
                break

        cycle_length = clock - old_clock
        regen_cycle_length_list.append(cycle_length)
        old_clock = clock
        regens += 1

        # TODO: result statistics
        # Insert the section that calculates the result statistics here

        job_A, job_B, job_C = [], [], []

    time_of_next_arrival = arrival_time_of_job(job_id)
    event = next_event(time_of_next_arrival)

    if event == "arrival":
        if curr_jobs > theta():
            # TODO: how is the arrival rate λ calculated exactly?
            backed_jobs += 1
            status("Job left")
            continue

        clock = time_of_next_arrival
        disk_visits[job_id] = 0
        c = job_category()
        add_job_to_station("CPU", job_id, c)
        (job_A if c == "A" else job_B if c == "B" else job_C).append(job_id)
        curr_jobs += 1
        job_id += 1
    else:
        station = event
        i = station_index(station)
        clock = STATION_start_time[i] + STATION_remaining_time[i]
        match station:
            case "CPU":
                vij_DISK_category = vij[station_index("DISK")][category_index(c)]
                disk_visit_probability = vij_DISK_category / (vij_DISK_category + 1)

                # CHECK
                target_station = "DISK" if U() < disk_visit_probability else "OUT"
                if target_station == "DISK":
                    disk_visits[STATION_current_job[i]] += 1
                add_job_to_station(
                    target_station, STATION_current_job[i], STATION_current_category[i]
                )
                decrease_CPU_remaining_time(STATION_remaining_time[i])
                STATION_queue[i].popleft()
                if not STATION_queue[i]:
                    set_current_counters("CPU", None, None, None, None)
                    STATION_LAST_EMPTY[i] = clock
                    empty_STATION[i] = True
                else:
                    job, category, remaining_time = STATION_queue[i][0]
                    set_current_counters("CPU", job, category, clock, remaining_time)
            case "DISK":
                add_job_to_station(
                    "CPU", STATION_current_job[i], STATION_current_category[i]
                )
                load_job_from_queue(station)
            case "OUT":
                curr_jobs -= 1
                load_job_from_queue(station)

    status(event)


print(f"Regen cycles = {regens}")
t = list(range(len(plot_list[0])))
plt.hist([x + y for (x, y) in zip(plot_list[0], plot_list[1])], bins="auto")
plt.title("Histogram of jobs on the system")
plt.xlabel("Number of jobs")
plt.ylabel("Events with this number of jobs")
plt.show()


### CHECK ARRIVAL RATES
print("Simulation, True, Number of events")
for key in arrivals_dict.keys():
    print(
        np.round(key, 1),
        float(np.round(sum(arrivals_dict[key]) / len(arrivals_dict[key]), 1)),
        len(arrivals_dict[key]),
        sep=", ",
    )


### OUTPUT RESULTS

# TODO: result statistics
# Ζητούμενα: λ, λj, R, Rj, Ui


R_A = round(float(average(R_A_different) / average(num_of_A_jobs_in_regen_cycle)), 3)
print(f"Average response time for Category A: {R_A}")
R_B = round(float(average(R_B_different) / average(num_of_B_jobs_in_regen_cycle)), 3)
print(f"Average response time for Category B: {R_B}")
R_C = round(float(average(R_C_different) / average(num_of_C_jobs_in_regen_cycle)), 3)
print(f"Average response time for Category C: {R_C}")
U_CPU = round(float(1 - average(U_CPU_different) / average(regen_cycle_length_list)), 3)
print(f"CPU utilization: {U_CPU}")
U_DISK = round(
    float(1 - average(U_DISK_different) / average(regen_cycle_length_list)), 3
)
print(f"DISK utilization: {U_DISK}")
U_OUT = round(float(1 - average(U_OUT_different) / average(regen_cycle_length_list)), 3)
print(f"OUT utilization: {U_OUT}")


### DISK 'VISITS'
the_list = []
for key in disk_visits:
    the_list.append(disk_visits[key])

N = 50
bins = np.linspace(0, N, N + 1)
answer = plt.hist(the_list, bins=bins)
plt.title("Disk Visits")
plt.ylabel("Number of jobs")
plt.xlabel("Number of Visits")
plt.show()

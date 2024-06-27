# Simulation program

import numpy as np
from scipy.stats import erlang
from collections import deque
import matplotlib.pyplot as plt

stations: tuple[str] = ("CPU", "DISK", "OUT")
categories: tuple[str] = ("A", "B", "C")
station_index: dict[str, int] = {}
category_index: dict[str, int] = {}
for i, station in enumerate(stations):
    station_index[station] = i
for j, category in enumerate(categories):
    category_index[category] = j

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

arrival_rate: float = 256776 / 45100  # λ: total arrivals per total time frame
average_arrival_time: float = 1 / arrival_rate

Dij: tuple[tuple[float]] = (
    (D_CPU_A, D_CPU_B, D_CPU_C),
    (D_DISK_A, D_DISK_B, D_DISK_C),
    (D_OUT_A, D_OUT_B, D_OUT_C),
)
vij: tuple[tuple[float]] = (
    (v_CPU_A, v_CPU_B, v_CPU_C),
    (v_DISK_A, v_DISK_B, v_DISK_C),
    (v_OUT_A, v_OUT_B, v_OUT_C),
)
Sij: tuple[list[float]] = ([], [], [])
for i in range(len(stations)):
    for j in range(len(categories)):
        Sij[i].append(Dij[i][j] / vij[i][j])

# Result variables
# X_k: list of numbers that correspond to the average number of each category
# X_k_per_cycle: tuple of list of numbers that correspond to the average number of a category of each cycle
lamda_j: list[float] = [None, None, None]
lamda_j_per_cycle: tuple[list[float]] = ([None], [None], [None])
R_j: list[float] = [None, None, None]
R_j_per_cycle: tuple[list[float]] = ([None], [None], [None])
U_i: list[float] = [None, None, None]
U_i_per_cycle: tuple[list[float]] = ([None], [None], [None])


def Erlang_4_service_time(category: str) -> float:
    # mean = k / λ => λ = k / mean
    k = 4
    lamda = k / Sij[station_index["CPU"]][category_index[category]]
    return erlang.rvs(a=k, scale=1 / lamda)


def Exponential_service_time(station: str, category: str) -> float:
    # mean = 1 / λ
    mean = Sij[station_index[station]][category_index[category]]
    return np.random.exponential(scale=mean)


U = lambda: float(np.random.uniform())


def job_category() -> str:
    """
    This function returns a category based on the number of
    arrivals of said category per total arrivals.

    Returns:
        str: A string indicating the category of the job that arrived.
    """
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


# sys_state_per_event has 3 lists, one for each station. Every time an event happens, the number of jobs
# using the station, either waiting in its queue or not, is appended on the station's list.
sys_state_per_event: tuple[list[int]] = ([], [], [])


def status() -> None:
    """
    This function is used to updated the number of jobs in the system
    every time a new event happens.
    """
    for i, station in enumerate(stations):
        if station == "CPU":
            sys_state_per_event[i].append(len(STATION_queue[i]))
        else:
            sys_state_per_event[i].append(
                0 if empty_STATION[i] else len(STATION_queue[i]) + 1
            )
    return


# STATION_queue has 3 fields, one for each station and it stores the corresponding
# station's queue. Each element has the following structure:
# (job_id, job_category, job_remaining_time)
STATION_queue = (deque(), deque(), deque())

STATION_current_job: list[int] = [None, None, None]
STATION_current_category: list[str] = [None, None, None]
STATION_start_time: list[float] = [None, None, None]
STATION_remaining_time: list[float] = [None, None, None]

STATION_EMPTY_TIME: list[float] = [0, 0, 0]  # re-initialized for each cycle
STATION_LAST_EMPTY: list[float] = [0, 0, 0]
empty_STATION: list[bool] = [True, True, True]

STATION_empty_time_per_cycle = ([None], [None], [None])


def set_current_counters(
    station: str, job: int, category: str, start_time: float, remaining_time: float
) -> None:
    """
    This function updates the variables of a station that correspond to the current job,
    current category, start time and remaining time.

    Args:
        station (str): The station of which the corresponding variables will change.
        job (int): Id of the job that is currently served by the station. If the station
            if the CPU, which uses a Processor sharing policy, this is the id of the job
            that requires the shortest time.
        category (str): The category of the job.
        start_time (float): The time that the station will start processing the job.
        remaining_time (float): The required service time by the job.
    """
    i = station_index[station]
    STATION_current_job[i] = job
    STATION_current_category[i] = category
    STATION_start_time[i] = start_time
    STATION_remaining_time[i] = remaining_time
    return


def next_event(arrival_time: float) -> str:
    """
    This function calculates the event that will happen sooner.

    Args:
        arrival_time (float): The time that the next arrival will happen.

    Returns:
        str: Returns either the station that will finish serving a job or "arrival" if a
            job will arrive in the system before a station finishes serving a job.
    """
    valid_sum = lambda x, y: x + y if x is not None and y is not None else None

    min_value = arrival_time
    result = "arrival"
    for i, station in enumerate(stations):
        STATION_finish_time = valid_sum(
            STATION_start_time[i], STATION_remaining_time[i]
        )
        if STATION_finish_time is not None and min_value > STATION_finish_time:
            min_value = STATION_finish_time
            result = station

    return result


# arrival_time(job_id: int) -> arrival_time_of_job: float
arrival_time_of_job = lambda job_id: job_id * average_arrival_time
# departure_time_of_job: dict[job_id, departure_time_of_job]
departure_time_of_job: dict[int, float] = {}
# disk_visits: dict[job_id, current_number_of_visits]
disk_visits: dict[int, float] = {}

clock: float = 0
job_id: int = 1
curr_jobs: int = 0  # re-initialized for each cycle (automatically)


def decrease_CPU_remaining_time(time_passed: float) -> None:
    """
    This function decreases the remaining time of all the jobs in CPU's queue by the given
    time.

    Args:
        time_passed (float): The time that the remaining times of the jobs will be reduced by.
    """
    for index in range(len(STATION_queue[station_index["CPU"]])):
        STATION_queue[station_index["CPU"]][index][2] -= time_passed
    return


def add_job_to_station(station: str, job_id: int, job_category: str) -> None:
    """
    This function adds a job to a station's queue.

    Args:
        station (str): The station-destination of the job.
        job_id (int): The job id of the job that will be added.
        job_category (str): The category of the job that will be added.
    """
    global clock
    i = station_index[station]
    service_time = float(
        Erlang_4_service_time(job_category)
        if station == "CPU"
        else Exponential_service_time(station, job_category)
    )

    def add_job_to_CPU_queue(
        job_id: int, job_category: str, service_time: float
    ) -> None:
        """
        This function adds a job to the correct position of the CPU's queue. All the jobs
        in the CPU's queue are sorted in increasing order of remaining service time.

        Args:
            job_id (int): The job's id that will be inserted into the queue.
            job_category (str): The job's category that will be inserted into the queue.
            service_time (flaot): The service time needed by the job.
        """
        insert_position = len(STATION_queue[station_index["CPU"]])
        for index, value in enumerate(STATION_queue[station_index["CPU"]]):
            if value[2] > service_time:
                insert_position = index
                break
        STATION_queue[station_index["CPU"]].insert(
            insert_position, [job_id, job_category, service_time]
        )
        return

    match station:
        case "CPU":
            if STATION_queue[i]:
                decrease_CPU_remaining_time(clock - STATION_start_time[i])
            add_job_to_CPU_queue(job_id, job_category, service_time)
            job, category, remaining_time = STATION_queue[i][0]
            set_current_counters(station, job, category, clock, remaining_time)
        case _:
            if empty_STATION[i]:
                set_current_counters(station, job_id, job_category, clock, service_time)
            else:
                STATION_queue[i].append([job_id, job_category, service_time])

    if empty_STATION[i]:
        STATION_EMPTY_TIME[i] += clock - STATION_LAST_EMPTY[i]
        STATION_LAST_EMPTY[i] = clock
        empty_STATION[i] = False

    return


def load_job_from_queue(station: str) -> None:
    """
    This function takes the next job in a station's queue and puts it in position for
    calculation. This station can only be a station that uses FIFO policy.

    Args:
        station (str): The station on which the load will happen.
    """
    global clock
    i = station_index[station]
    if not STATION_queue[i]:
        set_current_counters(station, None, None, None, None)
        empty_STATION[i] = True
        STATION_LAST_EMPTY[i] = clock
    else:
        job, category, remaining_time = STATION_queue[i].popleft()
        set_current_counters(station, job, category, clock, remaining_time)

    return


previous_regen_point = 0
jobs_per_category: tuple[list[int]] = ([], [], [])  # re-initialized for each cycle

cycle_index: int = 0
cycles_length: list[int] = [None]
Qj_cycle: tuple[list[int]] = ([None], [None], [None])
# Qt: list that contains the number of jobs in the system each time an event happens (clears after each cycles ends)
Qt: list[int] = []  # re-initialized for each cycle
# yi: interval of Qt divided by i-th cycle's length
yi: list[float] = [None]
check_condition = lambda ratio: ratio < 0.1
cycles_num_to_check: int = 20  # Trust degree: 95%

# average1() computes the average value of a list disregarding its 1st element
average1 = lambda lst: sum(lst[1:]) / (len(lst) - 1) if len(lst) > 1 else None
average = lambda lst: sum(lst) / len(lst) if lst else None
sum1 = lambda lst: sum(lst[1:]) if lst else None

while cycle_index < 1000:
    is_system_in_initial_state = curr_jobs == 0
    if is_system_in_initial_state and cycle_index != 0:
        cycles_length.append(clock - previous_regen_point)
        yi.append(sum(Qt) / cycles_length[cycle_index])
        previous_regen_point = clock

        # Result statistics: lamda_j, R_j, U_i
        for j in range(len(categories)):
            Qj_cycle[j].append(len(jobs_per_category[j]))
            lamda_j_per_cycle[j].append(
                Qj_cycle[j][cycle_index] / cycles_length[cycle_index]
            )

            R_j_temp = []
            for job in jobs_per_category[j]:
                R_j_temp.append(departure_time_of_job[job] - arrival_time_of_job(job))
            R_j_per_cycle[j].append(
                sum(R_j_temp) / len(R_j_temp) if jobs_per_category[j] else 0
            )

        for i in range(len(stations)):
            U_i_per_cycle[i].append(
                1 - STATION_EMPTY_TIME[i] / cycles_length[cycle_index]
            )
            STATION_empty_time_per_cycle[i].append(STATION_EMPTY_TIME[i])
            STATION_EMPTY_TIME[i] = 0

        if cycle_index % cycles_num_to_check == 0 and cycle_index != 0:
            y_bar = average1(yi)
            c_bar = average1(cycles_length)
            n = cycle_index

            sy_2 = 0
            sc_2 = 0
            syc = 0
            for index in range(1, cycle_index):
                sy_2 += (yi[index] - y_bar) ** 2
                sc_2 += (cycles_length[index] - c_bar) ** 2
                syc += (yi[index] - y_bar) * (cycles_length[index] - c_bar)

            sy_2 /= n - 1
            sc_2 /= n - 1
            syc /= n - 1

            R = y_bar / c_bar
            s_2 = sy_2 - 2 * R * syc + (R**2) * sc_2

            trust_degree = 100 - (100 / cycles_num_to_check)
            a = (100 - trust_degree) * 0.01
            z1_a_2 = 1 - a / 2
            s = float(np.sqrt(s_2))
            confidence_interval = z1_a_2 * s / (c_bar * float(np.sqrt(n)))

            if check_condition(confidence_interval / R):
                break

        for j in range(len(categories)):
            jobs_per_category[j].clear()
        Qt.clear()
        cycle_index += 1

    if cycle_index == 0:
        cycle_index = 1

    time_of_arrival: float = arrival_time_of_job(job_id)
    event = next_event(time_of_arrival)

    if event == "arrival":
        clock = time_of_arrival
        disk_visits[job_id] = 0
        c: str = job_category()
        add_job_to_station("CPU", job_id, c)
        jobs_per_category[category_index[c]].append(job_id)
        curr_jobs += 1
        job_id += 1
    else:
        station = event
        i = station_index[station]
        clock = STATION_start_time[i] + STATION_remaining_time[i]
        match station:
            case "CPU":
                vij_DISK_category = vij[station_index["DISK"]][category_index[c]]
                disk_visit_probability = vij_DISK_category / (vij_DISK_category + 1)

                target_station = "DISK" if U() < disk_visit_probability else "OUT"
                if target_station == "DISK":
                    disk_visits[STATION_current_job[i]] += 1
                add_job_to_station(
                    target_station, STATION_current_job[i], STATION_current_category[i]
                )
                decrease_CPU_remaining_time(STATION_remaining_time[i])
                STATION_queue[i].popleft()
                if not STATION_queue[i]:
                    set_current_counters(station, None, None, None, None)
                    STATION_LAST_EMPTY[i] = clock
                    empty_STATION[i] = True
                else:
                    job, category, remaining_time = STATION_queue[i][0]
                    set_current_counters(station, job, category, clock, remaining_time)
            case "DISK":
                add_job_to_station(
                    "CPU", STATION_current_job[i], STATION_current_category[i]
                )
                load_job_from_queue(station)
            case "OUT":
                departure_time_of_job[STATION_current_job[i]] = clock
                curr_jobs -= 1
                load_job_from_queue(station)

    Qt.append(curr_jobs)
    status()


print(f"Regenaration cycles = {cycle_index}")
total_num_of_events = len(sys_state_per_event[0])
CPU_jobs_per_cycle, DISK_jobs_per_cycle, OUT_jobs_per_cycle = sys_state_per_event
total_jobs_per_event = [
    cpu + disk + out
    for cpu, disk, out in zip(
        CPU_jobs_per_cycle, DISK_jobs_per_cycle, OUT_jobs_per_cycle
    )
]

# Output results: λ, λj, R, Rj, Ui, balking percentage

# Charts
## Histogram
plt.hist(total_jobs_per_event, bins="auto")
plt.title("Histogram of jobs on the system")
plt.xlabel("Number of jobs")
plt.ylabel("Events with this number of jobs")
plt.show()

t = list(range(total_num_of_events))

## Line chart
plt.figure(figsize=(10, 6))
plt.plot(t, total_jobs_per_event, marker="o", linestyle="-")

plt.title("Total Number of Jobs in the System Over Time")
plt.xlabel("Time (discrete events)")
plt.ylabel("Total Number of Jobs")
plt.grid(True)
plt.show()

## Bar chart
plt.figure(figsize=(10, 6))
plt.bar(t, total_jobs_per_event, color="skyblue")
plt.title("Total Number of Jobs in the System Over Time")
plt.xlabel("Time (discrete events)")
plt.ylabel("Total Number of Jobs")
plt.grid(True)
plt.show()

## Area chart
CPU_cumulative = CPU_jobs_per_cycle
DISK_cumulative = [cpu + disk for cpu, disk in zip(CPU_cumulative, DISK_jobs_per_cycle)]
OUT_cumulative = [disk + out for disk, out in zip(DISK_cumulative, OUT_jobs_per_cycle)]
plt.figure(figsize=(10, 6))
plt.fill_between(t, 0, CPU_cumulative, label="CPU", color="skyblue", alpha=0.5)
plt.fill_between(
    t,
    CPU_cumulative,
    DISK_cumulative,
    label="DISK",
    color="lightgreen",
    alpha=0.5,
)
plt.fill_between(
    t, DISK_cumulative, OUT_cumulative, label="OUT", color="salmon", alpha=0.5
)
plt.title("Total Number of Jobs in the System Over Time")
plt.xlabel("Time (discrete events)")
plt.ylabel("Total Number of Jobs")
plt.grid(True)
plt.show()


for j, category in enumerate(categories):
    lamda_j[j] = sum1(lamda_j_per_cycle[j]) / cycle_index
    print(
        f"Average arrival rate of category {category}: {float(np.round(lamda_j[j], 4))}"
    )
print(f"Total average arrival rate: {float(np.round(1 / average_arrival_time, 4))}")

for j, category in enumerate(categories):
    R_j_per_cycle_without_zeros = [elem for elem in R_j_per_cycle[j] if elem != 0]
    R_j[j] = round(
        float((average1(R_j_per_cycle_without_zeros) / average1(Qj_cycle[j]))), 3
    )
    print(f"Average response time for Category {category}: {float(np.round(R_j[j],3))}")
R = sum(R_j)
print(f"Average response time: {R}")

for i, station in enumerate(stations):
    # U_i[i] = round(float(1 - average1(U_i_per_cycle[i]) / average1(cycles_length)), 3)
    U_i[i] = (clock - sum1(STATION_empty_time_per_cycle[i])) / clock
    print(f"{station} utilization: {float(np.round(U_i[i], 3))}")


# times_of_N_jobs is a dictionary that maps the total number of jobs N in the system
# with the number of events that happend while N number of jobs where in the system
times_of_N_jobs: dict[int, int] = {}

for total_jobs in total_jobs_per_event:
    if not total_jobs in times_of_N_jobs.keys():
        times_of_N_jobs[total_jobs] = 0
    times_of_N_jobs[total_jobs] += 1

average_num_of_jobs = 0
for index, key in enumerate(times_of_N_jobs):
    average_num_of_jobs += key * times_of_N_jobs[key]
average_num_of_jobs /= total_num_of_events

print(
    f"Average number of jobs in the system (K): {float(np.round(average_num_of_jobs, 3))}"
)

entered_jobs: int = job_id - 1


def K() -> int:
    """
    This function returns the number of possible jobs in the system based on the
    probability found from 'total_jobs_per_event' list.

    Returns:
        int: The number of jobs in the system
    """

    U_temp = U()
    while U_temp == 1.0:
        U_temp = U()
    return total_jobs_per_event[int(U_temp * (total_num_of_events))]


balked_jobs: int = 0
theta = lambda: np.random.normal(loc=12, scale=3)

for event_index in range(total_num_of_events):
    if K() > theta():
        balked_jobs += 1

average_num_of_balked_jobs = balked_jobs / total_num_of_events

balking_percentage = float(
    np.round(
        average_num_of_balked_jobs / (entered_jobs + average_num_of_balked_jobs) * 100,
        7,
    )
)
print(f"Percentage of jobs that 'balked' from the system: {balking_percentage}%")

import numpy as np
import diptest 
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def _centroid(store, pids):               # среднее по группе пациентов, поверх уже усредненного представления пациента
    return store.sigs(pids).mean(0)         

def _medoid(store, pids):
    c = _centroid(store, pids)
    sig = store.sigs(pids)
    return pids[int(np.argmin(np.linalg.norm(sig - c, axis = 1)))]     # выбираем реальную клетку, самую близкую к центроиду

def _available(store, label, exclude_studies, used = ()):   
    m = store.meta
    mask = (m.label == label) & (~m.study.isin(exclude_studies)) & (~m.index.isin(used))
    return m.index[mask].tolist()             # 

def _fps_quantile(store, pids, K, q=0.8):
    if len(pids) <= K:
        return list(pids)
    sig = store.sigs(pids)
    chosen = [int(np.argmin(np.linalg.norm(sig - sig.mean(0), axis = 1)))]   # старт  медоид
    rest = [i for i in range(len(pids)) if i != chosen[0]]                   # позиции в pids, не сами id
    while rest and len(chosen) < K:
        d = np.array([min(np.linalg.norm(sig[i] - sig[c]) for c in chosen) for i in rest])
        pick = rest[int(np.argmin(np.abs(d - np.quantile(d, q))))]           # ближе к квантилю, не к максимуму
        chosen.append(pick)
        rest.remove(pick)
    return [pids[i] for i in chosen]

def _near_quantile(store, pids, anchor, K, q = 0.1):    
    if len(pids) <= K:
        return list(pids)
    sig = store.sigs(pids)
    d_anchor = np.linalg.norm(sig - anchor, axis = 1)     # каждый кандидат до якоря
    thr = max(np.quantile(d_anchor, q), np.partition(d_anchor, K - 1)[K - 1])  # q-ракушка ИЛИ K ближайших, что шире,  ракушка всегда >= K
    shell = [i for i in range(len(pids)) if d_anchor[i] <= thr]  # ближайшие, из которых потом отберём отдалённые друг от друга
    chosen = [shell[int(np.argmin(d_anchor[shell]))]]     # старт — самый близкий к якорю
    rest = [i for i in shell if i != chosen[0]]
    while rest and len(chosen) < K:
        dd = np.array([min(np.linalg.norm(sig[i] - sig[c]) for c in chosen) for i in rest])
        pick = rest[int(np.argmax(dd))]                   # внутри ракушка макс покрытие, стенка уже гарантирует близость
        chosen.append(pick)
        rest.remove(pick)
    return [pids[i] for i in chosen]

def _prep(store, cohort_pids, used):        # общий пролог ручек
    return set(used or ()), set(store.meta.loc[cohort_pids, "study"])

def _cohort_medoids(store, pids, m):        # m прототипов подмножества когорты: денойзенные якоря (компромисс centroid <-> per-case)
    if m <= 0 or not pids:
        return []
    if len(pids) <= m:
        return list(pids)
    lab = KMeans(m, n_init = 5, random_state = 0).fit_predict(store.sigs(pids))
    return [_medoid(store, [p for p, l in zip(pids, lab) if l == c]) for c in range(m)]  #type:ignore

def _nearest(store, pids, anchor):          # ближайший кандидат к якорю (pids уже без used)
    if not pids:
        return None
    return pids[int(np.argmin(np.linalg.norm(store.sigs(pids) - anchor, axis = 1)))]

def diverse_refs(store, cohort_pids, K_dis, K_nor, q = 0.8, used = None, case = 1, control = 0):
    """Коровы на пляже"""
    used, exclude = _prep(store, cohort_pids, used)
    refs = []
    if K_dis:
        refs += _fps_quantile(store, _available(store, case, exclude, used | set(refs)), K_dis, q)      # подготовка к мультиклассу
    if K_nor:
        refs += _fps_quantile(store, _available(store, control, exclude, used | set(refs)), K_nor, q)
    return refs

def matched_refs(store, cohort_pids, K_case, K_control, used = None, case = 1, control = 0):
    """Свиньи на травке"""
    used, exclude = _prep(store, cohort_pids, used)
    coh = store.meta.loc[cohort_pids]
    cases = coh.index[coh.label == case].tolist()
    controls = coh.index[coh.label == control].tolist()
    refs = []
    for anc in _cohort_medoids(store, cases, K_control):        # контроли к регионам случаев
        pick = _nearest(store, _available(store, control, exclude, used | set(refs)), store.sigs([anc])[0])
        if pick: refs.append(pick)
    for anc in _cohort_medoids(store, controls, K_case):        # случаи к регионам контролей
        pick = _nearest(store, _available(store, case, exclude, used | set(refs)), store.sigs([anc])[0])
        if pick: refs.append(pick)
    return refs
    
    
    
    
def _paired_studies(store, exclude_studies, case = 1, control = 0):
    m = store.meta
    out = {}
    for s, g in m[~m.study.isin(exclude_studies)].groupby("study"):
        dis = g.index[g.label == case].tolist()
        nor = g.index[g.label == control].tolist()
        if dis and nor:
            out[s] = (dis, nor)
    return out




"""
Функции, ждущие своего часа
"""
def diverse_refs_paired(store, cohort_pids, K, q = 0.8, seed = 0, case = 1, control = 0):
    """
    Корова на пляже
    """
    cohort_studies = set(store.meta.loc[cohort_pids, "study"])   # хэш-таблица с пациентами нашей когорты
    pstud = _paired_studies(store, cohort_studies, case, control)
    if not pstud:
        return []
    names = list(pstud)
    cents = {s: _centroid(store, dis + norm) for s, (dis, norm) in pstud.items()}  # распаковывает кортежи и ключи
    gc = np.mean([cents[s] for s in names], 0)  # усредняем по строкам (исследрваниям). Эта точка это центр масс корпуса
    
    start = min(names, key = lambda s: np.linalg.norm(cents[s] - gc))
    order, rest = [start], [s for s in names if s != start]
    while rest and len(order) < K:
        d = np.array([min(np.linalg.norm(cents[s] - cents[o]) for o in order) for s in rest])
        pick = rest[int(np.argmin(np.abs(d - np.quantile(d, q))))]
        order.append(pick)
        rest.remove(pick)

    refs = []
    for s in order:                    # уже ровно K (или меньше, если корпус мал)
        dis, nor = pstud[s]
        refs += [_medoid(store, dis), _medoid(store, nor)]
    return refs
        
def representatives(store, pids, n_min = 12, alpha = 0.05, max_k = 3):      # доказано полезная, но данных нужного размера у нас маловато что бы интегрировать сейчас
    if len(pids) < n_min:
        return [_medoid(store, pids)]
    X = store.sigs(pids) 
    pc1 = X @ np.linalg.svd(X - X.mean(0), full_matrices= False)[2][0]  # Сначала центрируем облако точек у нуля, получаем (U, S, Vt), U - новые координаты, S - дисперсия на каждой оси, Vt (k = min(n, d), d), full_matrices = False для экономии, только k, берем только Vt и из Vt первую компоненту
    _, pval = diptest.diptest(pc1)
    if pval >= alpha:              # проверка на унимодальность
        return [_medoid(store, pids)]
    best_k, best_s, best_lab = 1, -1, None
    for k in range(2, min(max_k, len(pids) - 1) + 1):
        lab = KMeans(k, n_init = 5, random_state=0).fit_predict(X)    # результаты кластеризации
        s = silhouette_score(X, lab)  # насколько каждая точка лучше вписывается в свой кластер, чем в другие, значение от -1 (ужасно) до 0(погранично) до 1 (хорошо)
        if s > best_s:
            best_k, best_s, best_lab = k, s, lab
    reps = []
    for c in range(best_k):
        members = [p for p, l in zip(pids, best_lab) if l == c]   #type:ignore разбираем паицентов по кластерам
        reps.append(_medoid(store, members))
    return reps 